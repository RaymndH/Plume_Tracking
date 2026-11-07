"""
Probabilistic smoke occupancy tracker — Task 16
Bayesian voxel-wise occupancy grid updated across 14 KRGX scans.
Centroid derived from probability expectation, not hard clustering.
"""
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime

# ── config ────────────────────────────────────────────────────────────────
OUT_DIR   = Path('krgx_hourly')
ALPHA     = 0.60    # shorter memory — column position updates faster,
                     # less accumulation of detached plume evidence
RHV_MAX   = 0.90    # smoke RhoHV threshold (matches krgx_hourly.py)
DBZ_MIN   = 10.0    # dBZ floor (matches DBZ_THRESHOLD)

# Likelihood values (probability that this voxel contains smoke)
P_SMOKE       = 0.92   # rhv < RHV_MAX AND dbz > DBZ_MIN
P_NON_SMOKE   = 0.08   # data present, thresholds not met
P_NO_DATA     = 0.50   # masked / no observation — neutral

# ── load scan list ────────────────────────────────────────────────────────
with open(OUT_DIR / 'hourly_results.json') as f:
    results = json.load(f)

grid_files = sorted(OUT_DIR.glob('grids_*.npz'))
print(f"Found {len(grid_files)} grid files, {len(results)} result entries")

# ── initialize occupancy grid from first scan's shape ────────────────────
d0       = np.load(grid_files[0])
nx       = len(d0['x_bins']) - 1
ny       = len(d0['y_bins']) - 1
nz       = len(d0['z_bins']) - 1
x_km     = (d0['x_bins'][:-1] + d0['x_bins'][1:]) / 2
y_km     = (d0['y_bins'][:-1] + d0['y_bins'][1:]) / 2
z_bins   = d0['z_bins']
z_km     = (z_bins[:-1] + z_bins[1:]) / 2
RADAR_ALT_KM = 2.530
z_km_msl = z_km + RADAR_ALT_KM

prob_grid = np.full((nx, ny, nz), 0.50, dtype=np.float32)  # uniform prior

# coordinate grids for expectation calculation
XX, YY, ZZ = np.meshgrid(x_km, y_km, z_km_msl, indexing='ij')

# ── Gaussian spatial prior centered on Bug Fire anchor ───────────────────
ANCHOR_X_KM  = -36.0
ANCHOR_Y_KM  =  -2.0
PRIOR_SIGMA_KM = 8.0   # e-folding radius — 8 km covers the column base
PRIOR_MAX    = 0.70    # prior probability AT the anchor
PRIOR_MIN    = 0.10    # prior probability far from anchor

# Build 2D Gaussian in x/y (same for all z levels)
dist2_xy = (XX - ANCHOR_X_KM)**2 + (YY - ANCHOR_Y_KM)**2   # (nx,ny,nz)
gauss_xy = np.exp(-dist2_xy / (2 * PRIOR_SIGMA_KM**2))      # (nx,ny,nz)

# Scale to [PRIOR_MIN, PRIOR_MAX]
prior_grid = (PRIOR_MAX - PRIOR_MIN) * gauss_xy + PRIOR_MIN
prior_grid = prior_grid.astype(np.float32)

_anchor_ix = np.argmin(np.abs(x_km - ANCHOR_X_KM))
_anchor_iy = np.argmin(np.abs(y_km - ANCHOR_Y_KM))
print(f"Prior grid: min={prior_grid.min():.3f}  "
      f"max={prior_grid.max():.3f}  "
      f"at anchor: {prior_grid[_anchor_ix, _anchor_iy, 0]:.3f}")

# ── output accumulators ───────────────────────────────────────────────────
times      = []
x_cents    = []
y_cents    = []
z_cents    = []
z_tops     = []   # 95th percentile of prob-weighted z
confidences = []  # total probability mass (proxy for detection strength)
scan_labels = []

# ── main loop ─────────────────────────────────────────────────────────────
for gf in grid_files:
    # parse scan time from filename: grids_YYYYMMDD_HHMMSS.npz
    stem = gf.stem  # e.g. grids_20260809_235957
    dt_str = stem.replace('grids_', '')
    dt = datetime.strptime(dt_str, '%Y%m%d_%H%M%S')

    d         = np.load(gf)
    rhv       = d['rhv_grid'].astype(float)
    dbz       = d['dbz_smooth'].astype(float)

    # ── build likelihood field ────────────────────────────────────────────
    has_data  = np.isfinite(rhv) & np.isfinite(dbz) & (dbz > 0)
    is_smoke  = has_data & (rhv < RHV_MAX) & (dbz > DBZ_MIN)
    is_clear  = has_data & ~is_smoke

    likelihood = np.full((nx, ny, nz), P_NO_DATA, dtype=np.float32)
    likelihood[is_smoke] = P_SMOKE
    likelihood[is_clear] = P_NON_SMOKE

    # ── Bayesian update with temporal decay ───────────────────────────────
    # 1. Decay prior toward neutral (0.5) — old evidence fades
    decay_target = np.where(has_data, prior_grid, 0.5)
    prob_grid = ALPHA * prob_grid + (1 - ALPHA) * decay_target

    # 2. Bayesian update: P_new ∝ P_prior * likelihood
    unnorm    = prob_grid * likelihood
    # Normalise per-voxel relative to complement
    # (binary: smoke vs non-smoke)
    complement = (1 - prob_grid) * (1 - likelihood)
    denom     = unnorm + complement
    prob_grid = np.where(denom > 0, unnorm / denom, 0.5).astype(np.float32)

    # ── derive centroid from probability expectation ──────────────────────
    # Use voxels above a soft threshold to avoid noise floor dominating
    weight    = np.maximum(prob_grid - 0.5, 0.0)   # only above-neutral mass
    total_w   = weight.sum()

    if total_w > 0:
        cx = float((XX * weight).sum() / total_w)
        cy = float((YY * weight).sum() / total_w)
        cz = float((ZZ * weight).sum() / total_w)

        # z_top: 95th percentile of probability-weighted z distribution
        # flatten and sort by z, accumulate weight
        z_flat = ZZ.ravel()
        w_flat = weight.ravel()
        sort_i = np.argsort(z_flat)
        cumw   = np.cumsum(w_flat[sort_i])
        cumw  /= cumw[-1]
        z_top  = float(z_flat[sort_i][np.searchsorted(cumw, 0.95)])
    else:
        cx = cy = cz = z_top = np.nan

    confidence = float(total_w)

    times.append(dt)
    x_cents.append(cx)
    y_cents.append(cy)
    z_cents.append(cz)
    z_tops.append(z_top)
    confidences.append(confidence)
    scan_labels.append(dt.strftime('%H:%M'))

    print(f"  {dt.strftime('%H:%MZ')}  "
          f"centroid=({cx:.1f}, {cy:.1f}, {cz:.2f} km MSL)  "
          f"z_top={z_top:.2f}  conf={confidence:.0f}")

# ── compute scan-to-scan drift ────────────────────────────────────────────
x_arr  = np.array(x_cents)
y_arr  = np.array(y_cents)
t_min  = np.array([(t - times[0]).total_seconds()/60 for t in times])

dx     = np.diff(x_arr)
dy     = np.diff(y_arr)
dt_min = np.diff(t_min)
speed  = np.sqrt(dx**2 + dy**2) / (dt_min / 60)
bearing = (np.degrees(np.arctan2(dx, dy)) % 360)

print(f"\nDrift statistics (all {len(speed)} intervals):")
print(f"  Speed   — mean={np.nanmean(speed):.1f}  "
      f"median={np.nanmedian(speed):.1f}  "
      f"max={np.nanmax(speed):.1f} km/hr")
print(f"  Bearing — mean={np.nanmean(bearing):.0f}°  "
      f"median={np.nanmedian(bearing):.0f}°")

# ── plots ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(12, 9))

# 1. Plan view centroid track
ax = axes[0, 0]
sc = ax.scatter(x_arr, y_arr, c=t_min, cmap='plasma', s=60, zorder=3)
ax.plot(x_arr, y_arr, 'k-', lw=1.2, zorder=2)
for i, lbl in enumerate(scan_labels):
    ax.annotate(lbl, (x_arr[i], y_arr[i]),
                fontsize=7, xytext=(4,4), textcoords='offset points')
ax.scatter(0, 0, marker='+', c='black', s=120, zorder=5, label='KRGX')
ax.scatter(-35.8, 1.1, marker='*', c='red', s=150, zorder=5,
           label='Bug Fire anchor')
plt.colorbar(sc, ax=ax, label='Minutes from first scan')
ax.set(xlabel='X km', ylabel='Y km',
       title='Probabilistic centroid track (plan view)', aspect='equal')
ax.legend(fontsize=8)

# 2. Confidence (probability mass) vs time
ax = axes[0, 1]
ax.plot(t_min, confidences, 'o-', color='steelblue')
ax.set(xlabel='Minutes from first scan',
       ylabel='Probability mass (above-neutral)',
       title='Detection confidence vs time')

# 3. Height vs time
ax = axes[1, 0]
ax.plot(t_min, z_cents, 'o-', color='darkorange', label='Z centroid')
ax.plot(t_min, z_tops,  's--', color='firebrick',  label='Z top (p95)')
ax.set(xlabel='Minutes from first scan', ylabel='km MSL',
       title='Plume height vs time (probabilistic)')
ax.legend(fontsize=8)

# 4. Scan-to-scan speed
ax = axes[1, 1]
t_mid = (t_min[:-1] + t_min[1:]) / 2
ax.bar(t_mid, speed, width=3, color='purple', alpha=0.7)
ax.axhline(20, color='red', lw=1.5, linestyle='--', label='20 km/hr ref')
ax.set(xlabel='Minutes from first scan', ylabel='km/hr',
       title='Scan-to-scan centroid speed')
ax.legend(fontsize=8)

plt.suptitle('KRGX — Bug Fire probabilistic plume tracker\n'
             '2026-08-09 23:59Z to 2026-08-10 00:56Z', fontsize=11)
plt.tight_layout()
out_png = OUT_DIR / 'task16_prob_tracker.png'
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print(f'\nSaved {out_png}')
