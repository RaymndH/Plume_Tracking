"""
KRGX hourly plume analysis
- Fetches one hour of NEXRAD Level-II scans
- For each scan: grids via histogramdd, smooths, finds hotspots via watershed
- Saves a 2D column-max PNG per scan with watershed regions overlaid
- Stores results list for further analysis
"""

import io, re, gzip, warnings, time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from itertools import combinations

import boto3
from botocore import UNSIGNED
from botocore.config import Config
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pyart
from skimage.morphology import skeletonize
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
from scipy.ndimage import gaussian_filter, binary_dilation

warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════
# CONFIG
# ═══════════════════════════════════════════════════════

SITE        = 'KRGX'
RADAR_LAT   = 39.7541
RADAR_LON   = -119.4608
RADAR_ALT_M = 2530.0

T_CENTER = datetime(2026, 8, 9, 23, 59, 57, tzinfo=timezone.utc)
T_BEFORE = timedelta(minutes=0)
T_AFTER  = timedelta(minutes=60)

X_MIN, X_MAX = -50_000, 10_000
Y_MIN, Y_MAX = -10_000,  25_000
Z_MAX        =   8_000

DXY, DZ = 250, 250

ELEV_MIN, ELEV_MAX = -0.5, 19.5

DBZ_THRESHOLD = 10
GAUSS_SIGMA   = 1.0

SADDLE_DROP_THRESH = 3.0   # dBZ — below this, merge two regions
PEAK_FRAC          = 0.65  # peak must be >= this fraction of global max
NEIGHBOR_RADIUS_M  = 2000  # m — neighborhood radius for peak validation
DIST_THRESHOLD_KM  = 20.0  # hotspots beyond this → unconfirmed source

OUT_DIR   = Path('./krgx_hourly')
OUT_DIR.mkdir(exist_ok=True)
CACHE_DIR = Path('./krgx_cache')
CACHE_DIR.mkdir(exist_ok=True)
BUCKETS   = ['unidata-nexrad-level2', 'noaa-nexrad-level2']

from pyproj import Proj, Transformer as _Transformer

_aeqd    = Proj(proj='aeqd', lat_0=RADAR_LAT, lon_0=RADAR_LON, units='m')
_wgs84   = Proj(proj='latlong', datum='WGS84')
_to_aeqd = _Transformer.from_proj(_wgs84, _aeqd, always_xy=True)

# GOES power-weighted centroid, Task 08b — independently verified anchor
FIRE_SOURCES = [
    (39.7631, -119.8788, 'Bug Fire (main)'),
]
fire_source_xy = []
for (_lat, _lon, _name) in FIRE_SOURCES:
    _xm, _ym = _to_aeqd.transform(_lon, _lat)
    fire_source_xy.append((_xm/1e3, _ym/1e3, _name))

# ═══════════════════════════════════════════════════════
# S3 HELPERS
# ═══════════════════════════════════════════════════════

s3     = boto3.client('s3', config=Config(signature_version=UNSIGNED),
                      region_name='us-east-1')
KEY_RE = re.compile(rf'{SITE}(\d{{8}})_(\d{{6}})_V06(\.gz)?$')

t_start = T_CENTER - T_BEFORE
t_end   = T_CENTER + T_AFTER


def list_keys(site, t_start, t_end):
    for bucket in BUCKETS:
        results, seen = [], set()
        for dt_anchor in [t_start, t_end]:
            day = dt_anchor.strftime('%Y/%m/%d')
            if day in seen: continue
            seen.add(day)
            try:
                pager = s3.get_paginator('list_objects_v2')
                for page in pager.paginate(Bucket=bucket, Prefix=f'{day}/{site}/'):
                    for obj in page.get('Contents', []):
                        k = obj['Key']
                        fname = k.split('/')[-1]
                        if '_MDM' in fname: continue
                        m = KEY_RE.match(fname)
                        if not m: continue
                        dt_f = datetime.strptime(
                            m.group(1) + m.group(2), '%Y%m%d%H%M%S'
                        ).replace(tzinfo=timezone.utc)
                        if t_start <= dt_f <= t_end:
                            results.append((bucket, k, dt_f))
            except: continue
        if results:
            results.sort(key=lambda x: x[2])
            print(f'  {len(results)} scans in {bucket}')
            return results
    return []


def fetch_raw(bucket, key):
    fname = key.split('/')[-1].removesuffix('.gz')
    cp = CACHE_DIR / fname
    if cp.exists(): return cp.read_bytes()
    for k in [key, key.removesuffix('.gz'), key.removesuffix('.gz') + '.gz']:
        try:
            raw = s3.get_object(Bucket=bucket, Key=k)['Body'].read()
            if k.endswith('.gz'): raw = gzip.decompress(raw)
            cp.write_bytes(raw)
            return raw
        except: continue
    raise FileNotFoundError(key)


# ═══════════════════════════════════════════════════════
# PER-SCAN PROCESSING
# ═══════════════════════════════════════════════════════

def process_scan(bucket, key, fname_dt):
    """
    Returns dict with hotspot info, or None if no echo found.
    """
    # ── fetch + read ──────────────────────────────────
    raw   = fetch_raw(bucket, key)
    radar = pyart.io.read_nexrad_archive(
        io.BytesIO(raw), include_fields=['reflectivity', 'cross_correlation_ratio']
    )
    keep = [
        i for i in range(radar.nsweeps)
        if ELEV_MIN <= float(radar.elevation['data'][radar.get_slice(i)].mean()) <= ELEV_MAX
    ]

    # ── extract gates ─────────────────────────────────
    gate_x, gate_y, gate_z, gate_v, gate_rhv = [], [], [], [], []
    for isw in keep:
        sl   = radar.get_slice(isw)
        refl = radar.fields['reflectivity']['data'][sl]
        if np.ma.is_masked(refl): refl = refl.filled(np.nan)
        refl   = refl.astype(float)
        gx     = radar.gate_x['data'][sl] / 1e3
        gy     = radar.gate_y['data'][sl] / 1e3
        gz_agl = radar.gate_altitude['data'][sl] / 1e3 - RADAR_ALT_M / 1e3
        valid  = (np.isfinite(refl) &
                  (gx     >= X_MIN/1e3) & (gx     <= X_MAX/1e3) &
                  (gy     >= Y_MIN/1e3) & (gy     <= Y_MAX/1e3) &
                  (gz_agl >= 0))
        gate_x.append(gx[valid].ravel())
        gate_y.append(gy[valid].ravel())
        gate_z.append(gz_agl[valid].ravel())
        gate_v.append(refl[valid].ravel())

        rhv_field = radar.fields.get('cross_correlation_ratio')
        if rhv_field is not None:
            rhv_data = rhv_field['data'][sl]
            rhv_arr  = np.where(
                np.ma.getmaskarray(rhv_data),
                np.nan,
                rhv_data.astype(float)
            )
            gate_rhv.append(rhv_arr[valid].ravel())
        else:
            gate_rhv.append(np.full(valid.sum(), np.nan))
    del radar, raw

    if not gate_x or sum(len(g) for g in gate_x) == 0:
        return None

    gx_all = np.concatenate(gate_x)
    gy_all = np.concatenate(gate_y)
    gz_all = np.concatenate(gate_z)
    gv_all = np.concatenate(gate_v)

    # ── histogramdd grid ──────────────────────────────
    x_bins = np.arange(X_MIN, X_MAX + DXY, DXY) / 1e3
    y_bins = np.arange(Y_MIN, Y_MAX + DXY, DXY) / 1e3
    z_bins = np.arange(0, Z_MAX/1e3 - RADAR_ALT_M/1e3 + DZ/1e3, DZ/1e3)

    pts        = np.column_stack([gx_all, gy_all, gz_all])
    H,      _  = np.histogramdd(pts, bins=[x_bins, y_bins, z_bins], weights=gv_all)
    counts, _  = np.histogramdd(pts, bins=[x_bins, y_bins, z_bins])
    with np.errstate(invalid='ignore'):
        dbz_grid = np.where(counts > 0, H / counts, np.nan)

    # RhoHV grid — per-field count to avoid bias from differential masking
    grhv_all   = np.concatenate(gate_rhv)
    rhv_finite = np.isfinite(grhv_all)
    if rhv_finite.sum() > 0:
        rhv_sum, _    = np.histogramdd(
            pts[rhv_finite], bins=[x_bins, y_bins, z_bins],
            weights=grhv_all[rhv_finite]
        )
        rhv_counts, _ = np.histogramdd(
            pts[rhv_finite], bins=[x_bins, y_bins, z_bins]
        )
        with np.errstate(invalid='ignore'):
            rhv_grid = np.where(rhv_counts > 0, rhv_sum / rhv_counts, np.nan)
    else:
        rhv_grid = np.full_like(dbz_grid, np.nan)

    x_km     = (x_bins[:-1] + x_bins[1:]) / 2
    y_km     = (y_bins[:-1] + y_bins[1:]) / 2
    z_km     = (z_bins[:-1] + z_bins[1:]) / 2
    z_km_msl = z_km + RADAR_ALT_M / 1e3

    # ── smooth + binary ───────────────────────────────
    dbz_smooth = gaussian_filter(np.nan_to_num(dbz_grid, nan=0.0), sigma=GAUSS_SIGMA)
    binary_3d  = dbz_smooth > DBZ_THRESHOLD

    RHV_MAX       = 0.90
    smoke_mask_3d = (
        (dbz_smooth > DBZ_THRESHOLD) &
        np.isfinite(rhv_grid) &
        (rhv_grid < RHV_MAX)
    )
    smoke_vox = int(smoke_mask_3d.sum())
    smoke_footprint_2d = smoke_mask_3d.any(axis=2)   # (nx, ny) bool

    if binary_3d.sum() < 10:
        return None

    # ── peak detection ────────────────────────────────
    col_max_2d       = dbz_smooth.max(axis=2)
    global_max       = dbz_smooth.max()
    MIN_PEAK_VAL     = global_max * PEAK_FRAC
    MIN_NEIGHBOR_VOX = round((NEIGHBOR_RADIUS_M / DXY) ** 2 * 0.3)

    coords_2d  = peak_local_max(col_max_2d, min_distance=5, threshold_abs=5)
    real_peaks = []
    threshold_abs = 5
    for px, py in coords_2d:
        peak_val = col_max_2d[px, py]
        r_vox    = round(NEIGHBOR_RADIUS_M / DXY)
        x0, x1  = max(0, px-r_vox), min(col_max_2d.shape[0], px+r_vox)
        y0, y1  = max(0, py-r_vox), min(col_max_2d.shape[1], py+r_vox)
        n_above  = (col_max_2d[x0:x1, y0:y1] > DBZ_THRESHOLD).sum()
        if n_above >= MIN_NEIGHBOR_VOX and peak_val >= threshold_abs:
            real_peaks.append((px, py, peak_val))

    real_peaks.sort(key=lambda x: -x[2])
    chosen = [(px, py) for px, py, _ in real_peaks]

    if not chosen:
        return None

    # ── watershed ─────────────────────────────────────
    markers = np.zeros_like(dbz_smooth, dtype=int)
    for i, (px, py) in enumerate(chosen):
        pz = np.argmax(dbz_smooth[px, py, :])
        markers[px, py, pz] = i + 1

    labels_ws = watershed(-dbz_smooth, markers, mask=binary_3d)

    # ── saddle-based merging ──────────────────────────
    merge_map = {i+1: i+1 for i in range(len(chosen))}
    for r1, r2 in combinations(range(1, len(chosen)+1), 2):
        reg1   = labels_ws == r1
        reg2   = labels_ws == r2
        border = binary_dilation(reg1) & binary_dilation(reg2) & binary_3d
        if border.sum() == 0: continue
        saddle   = dbz_smooth[border].max()
        peak1    = dbz_smooth[reg1].max()
        peak2    = dbz_smooth[reg2].max()
        drop_min = min(peak1 - saddle, peak2 - saddle)
        if drop_min <= SADDLE_DROP_THRESH:
            old = merge_map[r2]
            new = merge_map[r1]
            for k in merge_map:
                if merge_map[k] == old:
                    merge_map[k] = new

    groups = {}
    for rid, gid in merge_map.items():
        region = labels_ws == rid
        groups.setdefault(gid, np.zeros_like(binary_3d, dtype=bool))
        groups[gid] |= region

    # ── collect hotspot stats ─────────────────────────
    hotspots = []
    for gid, region in groups.items():
        if region.sum() == 0: continue
        ix, iy, iz = np.where(region)
        cx, cy = float(x_km[ix].mean()), float(y_km[iy].mean())
        dists  = [np.sqrt((cx-sx)**2 + (cy-sy)**2)
                  for (sx, sy, _) in fire_source_xy]
        ni     = int(np.argmin(dists))
        ndist  = round(dists[ni], 2)
        if ndist <= DIST_THRESHOLD_KM:
            src_name      = fire_source_xy[ni][2]
            src_confirmed = True
        else:
            src_name      = f'Unconfirmed (nearest: {fire_source_xy[ni][2]} {ndist:.1f} km)'
            src_confirmed = False
        hotspots.append({
            'x_mean':  cx,
            'y_mean':  cy,
            'z_mean':  float(z_km_msl[iz].mean()),
            'z_max':   float(z_km_msl[iz].max()),
            'dbz_max': float(dbz_smooth[region].max()),
            'voxels':  int(region.sum()),
            'region':  region,
            'source_name':      src_name,
            'source_dist_km':   ndist,
            'source_confirmed': src_confirmed,
            'smoke_vox':        int((smoke_mask_3d & region).sum()),
            'smoke_frac':       round(
                (smoke_mask_3d & region).sum() / max(region.sum(), 1), 3),
        })
    hotspots.sort(key=lambda h: -h['dbz_max'])

    # ── 2D plot ───────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 7))
    norm = mcolors.Normalize(vmin=0, vmax=25)
    ax.pcolormesh(x_km, y_km, col_max_2d.T, cmap='jet', norm=norm)
    plt.colorbar(ax.collections[0], ax=ax, label='dBZ')

    colors = plt.cm.tab10(np.linspace(0, 1, max(len(hotspots), 1)))
    for ci, hs in enumerate(hotspots):
        region    = hs['region']
        footprint = region.any(axis=2).T
        ax.contour(x_km, y_km, footprint.astype(float),
                   levels=[0.5], colors=[colors[ci]], linewidths=2)
        ax.scatter(hs['x_mean'], hs['y_mean'], c=[colors[ci]],
                   s=100, marker='x', linewidths=2, zorder=5)
        ax.annotate(f'{hs["dbz_max"]:.1f} dBZ\n{hs["z_max"]:.2f} km top',
                    (hs['x_mean'], hs['y_mean']),
                    color=colors[ci], fontsize=8,
                    xytext=(5, 5), textcoords='offset points')

    smoke_footprint = smoke_mask_3d.any(axis=2).T
    if smoke_footprint.any():
        ax.contourf(x_km, y_km, smoke_footprint.astype(float),
                    levels=[0.5, 1.5], colors=['purple'], alpha=0.25,
                    zorder=2)
        ax.contour(x_km, y_km, smoke_footprint.astype(float),
                   levels=[0.5], colors=['purple'], linewidths=1,
                   linestyles='--', zorder=3)

    ax.scatter(0, 0, marker='+', c='white', s=120, zorder=6, label='KRGX')
    ax.set(aspect='equal', xlabel='X km', ylabel='Y km',
           title=f'KRGX  {fname_dt.strftime("%H:%M:%S Z")}  '
                 f'| {len(hotspots)} hotspot(s)')
    ax.legend(fontsize=8)
    plt.tight_layout()
    out_path = OUT_DIR / f'krgx_{fname_dt.strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()

    return {
        'dt':               fname_dt,
        'hotspots':         [{k: v for k, v in h.items() if k != 'region'}
                              for h in hotspots],
        'png':              out_path,
        'smoke_vox_total':  smoke_vox,
        'smoke_footprint_2d': smoke_footprint_2d,
    }


# ═══════════════════════════════════════════════════════
# MAIN LOOP
# ═══════════════════════════════════════════════════════

print('Listing scans...')
scan_list = list_keys(SITE, t_start, t_end)
for _, key, dt in scan_list:
    print(f'  {dt.strftime("%H:%M:%SZ")}  {key.split("/")[-1]}')

results = []

for si, (bucket, key, fname_dt) in enumerate(scan_list):
    t0 = time.perf_counter()
    print(f'\n[{si+1}/{len(scan_list)}] {fname_dt.strftime("%H:%M:%SZ")}')
    try:
        res = process_scan(bucket, key, fname_dt)
        if res is None:
            print('  no echo / insufficient coverage')
            continue
        for hi, hs in enumerate(res['hotspots']):
            print(f'  Hotspot {hi+1}: x={hs["x_mean"]:.1f} y={hs["y_mean"]:.1f} '
                  f'z_mean={hs["z_mean"]:.2f} z_top={hs["z_max"]:.2f} km MSL  '
                  f'max={hs["dbz_max"]:.1f} dBZ  vox={hs["voxels"]}')
        print(f'  Smoke mask: {res["smoke_vox_total"]} vox total')
        for hi, hs in enumerate(res['hotspots']):
            confirmed_str = '' if hs['source_confirmed'] else ' [UNCONFIRMED]'
            print(f'    Hotspot {hi+1} -> {hs["source_name"]}'
                  f'{confirmed_str}  '
                  f'dist={hs["source_dist_km"]:.1f} km  '
                  f'smoke={hs["smoke_frac"]*100:.0f}%')
        print(f'  saved {res["png"].name}  ({time.perf_counter()-t0:.1f}s)')
        results.append(res)
    except Exception as e:
        print(f'  ERROR: {e}')

print(f'\nDone: {len(results)}/{len(scan_list)} scans processed')
print(f'PNGs saved to {OUT_DIR.resolve()}')

import json

def _serial(obj):
    """JSON serializer for numpy types."""
    if isinstance(obj, (np.integer,)): return int(obj)
    if isinstance(obj, (np.floating,)): return float(obj)
    if isinstance(obj, (np.ndarray,)): return obj.tolist()
    raise TypeError(f"Not serializable: {type(obj)}")

out_json = OUT_DIR / 'hourly_results.json'
with open(out_json, 'w') as f:
    json.dump([
        {
            'dt': r['dt'].isoformat(),
            'smoke_vox_total': r['smoke_vox_total'],
            'hotspots': r['hotspots'],   # already dicts, region excluded
        }
        for r in results
    ], f, indent=2, default=_serial)
print(f'Results saved to {out_json}')

# Save per-scan smoke footprints as compressed npz
footprints = {
    r['dt'].strftime('%Y%m%d_%H%M%S'): r['smoke_footprint_2d']
    for r in results
    if 'smoke_footprint_2d' in r
}
np.savez_compressed(OUT_DIR / 'smoke_footprints.npz', **footprints)
print(f'Footprints saved to {OUT_DIR}/smoke_footprints.npz')
