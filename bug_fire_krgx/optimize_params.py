"""
optimize_params.py
──────────────────
Automatically optimize key firelines algorithm parameters for a fire
by maximizing IoU against VIIRS perimeters.

Usage (from a fire run directory):
    python ..\..\optimize_params.py
    python ..\..\optimize_params.py --method grid      # default
    python ..\..\optimize_params.py --method bayesian  # requires optuna
    python ..\..\optimize_params.py --quick            # 2-param grid only

Output:
    best_params.yaml  — best parameter set found
    opt_results.csv   — all trials with scores
"""

import argparse
import copy
import csv
import glob
import os
import pickle
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import numpy as np
import yaml
from pyproj import Proj
from shapely.geometry import Polygon
from shapely.ops import unary_union

SCRIPT_DIR = Path(r"C:\Users\raymondharding\Research\Runs")


# ── helpers ───────────────────────────────────────────────────────────────────

def is_nil(v):
    return v is None or (isinstance(v, str) and v.strip().lower() == "nil")


def load_yaml(path):
    with open(path) as f:
        return yaml.safe_load(f)


def save_yaml(cfg, path):
    with open(path, "w") as f:
        yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)


def iou(poly_a, poly_b):
    """Intersection over Union of two shapely Polygons."""
    if poly_a is None or poly_b is None:
        return 0.0
    try:
        inter = poly_a.intersection(poly_b).area
        union = poly_a.union(poly_b).area
        return float(inter / union) if union > 0 else 0.0
    except Exception:
        return 0.0


def load_radar_perimeter_at(fire_cloud, target_dt, max_gap_s=1800):
    """
    Find the fire_cloud frame closest to target_dt and return its perimeter
    as a Shapely Polygon in radar-relative aeqd metres.
    """
    best = None
    best_gap = float("inf")
    for fc in fire_cloud:
        gap = abs((fc["dt"] - target_dt).total_seconds())
        if gap < best_gap:
            best_gap = gap
            best = fc
    if best is None or best_gap > max_gap_s:
        return None
    try:
        return Polygon(zip(best["xperim"], best["yperim"]))
    except Exception:
        return None


def load_viirs_perimeters(perimdir, proj):
    """
    Load perimeters.geojson, reproject to radar-relative aeqd,
    return list of (datetime, Polygon).
    """
    path = Path(perimdir) / "perimeters.geojson"
    if not path.exists():
        return []
    gdf = gpd.read_file(path)
    results = []
    for _, row in gdf.iterrows():
        dt_str = row.get("datetime")
        if not dt_str:
            continue
        try:
            dt = datetime.fromisoformat(str(dt_str))
            if dt.tzinfo is None:
                dt = dt.replace(tzinfo=timezone.utc)
        except Exception:
            continue
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue
        if geom.geom_type == "MultiPolygon":
            geom = max(geom.geoms, key=lambda g: g.area)
        # reproject from WGS84 lon/lat to radar aeqd
        try:
            coords = np.array(geom.exterior.coords)
            xy = np.array([proj(lon, lat) for lon, lat in coords])
            poly = Polygon(xy)
            if poly.is_valid and not poly.is_empty:
                results.append((dt, poly))
        except Exception:
            continue
    return sorted(results, key=lambda x: x[0])


def score_trial(pkl_path, viirs_perims):
    """
    Compute mean IoU between radar perimeters and nearest VIIRS perimeters.
    """
    if not Path(pkl_path).exists():
        return 0.0
    try:
        with open(pkl_path, "rb") as f:
            fire_cloud = pickle.load(f)
    except Exception:
        return 0.0

    if not fire_cloud or not viirs_perims:
        # fall back to coverage score: fraction of frames with any detection
        n_det = sum(1 for fc in fire_cloud if len(fc.get("xkeep", [])) > 0)
        return float(n_det) / max(len(fire_cloud), 1)

    ious = []
    for vt, vpoly in viirs_perims:
        rpoly = load_radar_perimeter_at(fire_cloud, vt)
        ious.append(iou(rpoly, vpoly))

    # also penalize frames with zero detections
    n_det = sum(1 for fc in fire_cloud if len(fc.get("xkeep", [])) > 0)
    coverage = float(n_det) / max(len(fire_cloud), 1)

    mean_iou = float(np.mean(ious)) if ious else 0.0
    return 0.7 * mean_iou + 0.3 * coverage


def run_trial(cfg, trial_params, config_path, run_name):
    """
    Apply trial_params to cfg, write temp YAML, run radar_perimeters.py,
    return score.
    """
    trial_cfg = copy.deepcopy(cfg)
    r = trial_cfg["run"]
    for k, v in trial_params.items():
        r[k] = v
    r["run_name"] = f"opt_{run_name}"
    r["save_frames"] = False

    tmp_yaml = Path(config_path).parent / "_opt_trial.yaml"
    save_yaml(trial_cfg, tmp_yaml)

    t0 = time.perf_counter()
    result = subprocess.run(
        [sys.executable, str(SCRIPT_DIR / "radar_perimeters.py"), str(tmp_yaml)],
        cwd=Path(config_path).parent,
        capture_output=True, text=True, timeout=600
    )
    elapsed = time.perf_counter() - t0

    if result.returncode != 0:
        print(f"    [trial failed in {elapsed:.0f}s]")
        print(result.stderr[-500:] if result.stderr else "no stderr")
        return 0.0, elapsed

    pkl = Path(config_path).parent / r["output_dir"] / f"opt_{run_name}" / "fire_cloud_final.pkl"
    score = score_trial(pkl, viirs_perims)
    return score, elapsed


# ── search spaces ─────────────────────────────────────────────────────────────

GRID_QUICK = {
    "z75_floor":      [20, 23, 25, 27, 30],
    "outlier_thresh": [5, 10, 15, 20, 25],
}

GRID_FULL = {
    "z75_floor":      [20, 23, 25, 27, 30],
    "outlier_thresh": [5, 10, 15, 20],
    "dist_cutoff":    [2000, 3000, 5000],
    "outlier_weight": [0.0, 0.25, 0.5],
}


def grid_search(cfg, config_path, run_name, space, viirs_perims_ref):
    """Exhaustive grid search over parameter space."""
    global viirs_perims
    viirs_perims = viirs_perims_ref

    keys   = list(space.keys())
    values = list(space.values())
    from itertools import product
    combos = list(product(*values))
    print(f"Grid search: {len(combos)} trials over {keys}")

    results = []
    best_score = -1
    best_params = {}

    for i, combo in enumerate(combos):
        params = dict(zip(keys, combo))
        print(f"  [{i+1}/{len(combos)}] {params}", end="  ")
        score, elapsed = run_trial(cfg, params, config_path, run_name)
        print(f"→ score={score:.4f}  ({elapsed:.0f}s)")
        results.append({**params, "score": score, "elapsed": elapsed})
        if score > best_score:
            best_score = score
            best_params = params

    return best_params, best_score, results


def bayesian_search(cfg, config_path, run_name, viirs_perims_ref, n_trials=40):
    """Bayesian optimization using optuna."""
    try:
        import optuna
        optuna.logging.set_verbosity(optuna.logging.WARNING)
    except ImportError:
        print("optuna not installed — pip install optuna")
        sys.exit(1)

    global viirs_perims
    viirs_perims = viirs_perims_ref

    def objective(trial):
        params = {
            "z75_floor":      trial.suggest_int("z75_floor", 18, 32),
            "outlier_thresh": trial.suggest_int("outlier_thresh", 3, 30),
            "dist_cutoff":    trial.suggest_int("dist_cutoff", 1000, 8000, step=500),
            "outlier_weight": trial.suggest_float("outlier_weight", 0.0, 1.0, step=0.1),
            "sw_min_thresh":  trial.suggest_float("sw_min_thresh", 2.0, 5.0, step=0.5),
            "sw_max_thresh":  trial.suggest_float("sw_max_thresh", 4.0, 8.0, step=0.5),
        }
        score, _ = run_trial(cfg, params, config_path, run_name)
        return score

    study = optuna.create_study(direction="maximize",
                                 sampler=optuna.samplers.TPESampler(seed=42))
    study.optimize(objective, n_trials=n_trials, show_progress_bar=True)

    best = study.best_params
    results = [
        {**t.params, "score": t.value, "elapsed": 0}
        for t in study.trials
    ]
    return best, study.best_value, results


# ── main ──────────────────────────────────────────────────────────────────────

viirs_perims = []   # set by search functions

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", nargs="?", default=None)
    parser.add_argument("--method", choices=["grid", "bayesian"], default="grid")
    parser.add_argument("--quick",  action="store_true",
                        help="2-param grid only (fastest)")
    parser.add_argument("--trials", type=int, default=40,
                        help="number of Bayesian trials")
    args = parser.parse_args()

    config_path = args.config or glob.glob("*.yaml")[0]
    print(f"config: {config_path}")
    cfg = load_yaml(config_path)

    run_name = cfg.get("name", Path(config_path).stem)

    # build radar projection for VIIRS reprojection
    rs       = cfg["radar_station"]
    proj     = Proj(proj="aeqd", lat_0=float(rs["lat"]),
                    lon_0=float(rs["lon"]), datum="WGS84")
    perimdir = cfg.get("viirs", {}).get("perimdir", "./perims")
    global viirs_perims
    viirs_perims = load_viirs_perimeters(perimdir, proj)
    print(f"VIIRS perimeters: {len(viirs_perims)}")
    for vt, vp in viirs_perims:
        print(f"  {vt:%Y-%m-%d %H:%M}  area={vp.area/1e6:.1f} km²")

    # run optimization
    if args.method == "bayesian":
        best_params, best_score, results = bayesian_search(
            cfg, config_path, run_name, viirs_perims, n_trials=args.trials)
    else:
        space = GRID_QUICK if args.quick else GRID_FULL
        best_params, best_score, results = grid_search(
            cfg, config_path, run_name, space, viirs_perims)

    # report
    print(f"\n{'='*60}")
    print(f"Best score: {best_score:.4f}")
    print(f"Best params: {best_params}")

    # write best params back to a new YAML
    best_cfg = copy.deepcopy(cfg)
    for k, v in best_params.items():
        best_cfg["run"][k] = v
    best_cfg["run"]["run_name"] = None
    best_cfg["run"]["save_frames"] = False
    out_yaml = Path(config_path).parent / "best_params.yaml"
    save_yaml(best_cfg, out_yaml)
    print(f"Best params written → {out_yaml}")

    # write results CSV
    out_csv = Path(config_path).parent / "opt_results.csv"
    if results:
        with open(out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
            w.writeheader()
            w.writerows(results)
        print(f"All results → {out_csv}")

    # clean up temp yaml
    tmp = Path(config_path).parent / "_opt_trial.yaml"
    if tmp.exists():
        tmp.unlink()


if __name__ == "__main__":
    main()
