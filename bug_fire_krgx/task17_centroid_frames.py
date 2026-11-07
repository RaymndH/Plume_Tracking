"""
Per-timestep radar images with probabilistic centroid overlay — Task 17
Loads grids_*.npz (raw dBZ/RhoHV per scan) and task16 centroid output,
overlays centroid position on 2D column-max dBZ for each scan.
"""
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
from datetime import datetime

OUT_DIR   = Path('krgx_hourly')
FRAME_DIR = OUT_DIR / 'centroid_frames'
FRAME_DIR.mkdir(exist_ok=True)

RADAR_ALT_KM = 2.530
DBZ_MIN      = 10.0
RHV_MAX      = 0.90

# ── re-run the probabilistic tracker inline to get centroids ──────────────
# (copy the core loop from task16_prob_tracker.py — no import needed)
ALPHA       = 0.60    # shorter memory — column position updates faster,
                       # less accumulation of detached plume evidence
P_SMOKE     = 0.92
P_NON_SMOKE = 0.08
P_NO_DATA   = 0.50

grid_files = sorted(OUT_DIR.glob('grids_*.npz'))
d0         = np.load(grid_files[0])
nx  = len(d0['x_bins']) - 1
ny  = len(d0['y_bins']) - 1
nz  = len(d0['z_bins']) - 1
x_km     = (d0['x_bins'][:-1] + d0['x_bins'][1:]) / 2
y_km     = (d0['y_bins'][:-1] + d0['y_bins'][1:]) / 2
z_km_msl = (d0['z_bins'][:-1] + d0['z_bins'][1:]) / 2 + RADAR_ALT_KM

prob_grid = np.full((nx, ny, nz), 0.50, dtype=np.float32)
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

centroids   = []   # (dt, cx, cy, cz, confidence)

for gf in grid_files:
    stem  = gf.stem.replace('grids_', '')
    dt    = datetime.strptime(stem, '%Y%m%d_%H%M%S')
    d     = np.load(gf)
    rhv   = d['rhv_grid'].astype(float)
    dbz   = d['dbz_smooth'].astype(float)

    has_data  = np.isfinite(rhv) & np.isfinite(dbz) & (dbz > 0)
    is_smoke  = has_data & (rhv < RHV_MAX) & (dbz > DBZ_MIN)
    is_clear  = has_data & ~is_smoke

    likelihood = np.full((nx, ny, nz), P_NO_DATA, dtype=np.float32)
    likelihood[is_smoke] = P_SMOKE
    likelihood[is_clear] = P_NON_SMOKE

    decay_target = np.where(has_data, prior_grid, 0.5)
    prob_grid  = ALPHA * prob_grid + (1 - ALPHA) * decay_target
    unnorm     = prob_grid * likelihood
    complement = (1 - prob_grid) * (1 - likelihood)
    denom      = unnorm + complement
    prob_grid  = np.where(denom > 0, unnorm / denom, 0.5).astype(np.float32)

    weight  = np.maximum(prob_grid - 0.5, 0.0)
    total_w = weight.sum()

    if total_w > 0:
        cx = float((XX * weight).sum() / total_w)
        cy = float((YY * weight).sum() / total_w)
        cz = float((ZZ * weight).sum() / total_w)
    else:
        cx = cy = cz = np.nan

    centroids.append((dt, cx, cy, cz, float(total_w), dbz, rhv))

# ── per-scan frames ───────────────────────────────────────────────────────
norm     = mcolors.Normalize(vmin=0, vmax=35)
cmap_dbz = plt.cm.jet
smoke_cmap = mcolors.ListedColormap(['none', 'purple'])

# accumulate centroid trail for overlay
trail_x, trail_y = [], []

for i, (dt, cx, cy, cz, conf, dbz, rhv) in enumerate(centroids):
    trail_x.append(cx)
    trail_y.append(cy)

    # 2D column-max dBZ
    col_max = dbz.max(axis=2)   # (nx, ny)

    # 2D smoke mask footprint (column-any of prob > 0.6)
    smoke_fp = (prob_grid > 0.6).any(axis=2)   # approximate — use
    # the current prob_grid state isn't accessible here after the loop;
    # use rhv/dbz directly for the footprint overlay instead
    smoke_2d = (
        (np.isfinite(rhv) & (rhv < RHV_MAX) & (dbz > DBZ_MIN))
        .any(axis=2)
    )

    fig, ax = plt.subplots(figsize=(8, 7))

    # dBZ background
    pcm = ax.pcolormesh(x_km, y_km, col_max.T,
                        cmap=cmap_dbz, norm=norm, shading='auto')
    plt.colorbar(pcm, ax=ax, label='dBZ (column max)', pad=0.02)

    # smoke mask footprint (purple tint)
    ax.contourf(x_km, y_km, smoke_2d.T.astype(float),
                levels=[0.5, 1.5], colors=['purple'], alpha=0.20, zorder=2)
    ax.contour(x_km, y_km, smoke_2d.T.astype(float),
               levels=[0.5], colors=['purple'], linewidths=1,
               linestyles='--', zorder=3)

    # Prior boundary (where prior = 0.5 — tracker's neutral zone)
    prior_2d = prior_grid[:, :, 0]   # same at all z levels
    ax.contour(x_km, y_km, prior_2d.T,
               levels=[0.50], colors=['white'],
               linewidths=1.0, linestyles=':', alpha=0.6, zorder=3)

    # centroid trail (previous scans, faded)
    if len(trail_x) > 1:
        ax.plot(trail_x[:-1], trail_y[:-1],
                'w--', lw=1.0, alpha=0.5, zorder=4)
        ax.scatter(trail_x[:-1], trail_y[:-1],
                   c='white', s=20, alpha=0.4, zorder=4)

    # current centroid
    if np.isfinite(cx):
        ax.scatter(cx, cy, c='yellow', s=180, marker='*',
                   edgecolors='black', linewidths=0.8,
                   zorder=6, label=f'Centroid z={cz:.2f} km MSL')
        ax.annotate(
            f'{dt.strftime("%H:%MZ")}\n({cx:.1f}, {cy:.1f}) km\n'
            f'z={cz:.2f} km MSL\nconf={conf:.0f}',
            xy=(cx, cy),
            xytext=(8, 8), textcoords='offset points',
            fontsize=8, color='yellow',
            bbox=dict(boxstyle='round,pad=0.2', fc='black', alpha=0.5),
            zorder=7
        )

    # Bug Fire anchor
    ax.scatter(-36.0, -2.0, marker='*', c='red', s=150,
               edgecolors='white', linewidths=0.5,
               zorder=5, label='Bug Fire anchor')

    # radar site
    ax.scatter(0, 0, marker='+', c='white', s=120, zorder=5, label='KRGX')

    ax.set(aspect='equal',
           xlim=(x_km.min(), x_km.max()),
           ylim=(y_km.min(), y_km.max()),
           xlabel='X km (aeqd from KRGX)',
           ylabel='Y km (aeqd from KRGX)',
           title=f'KRGX  {dt.strftime("%Y-%m-%d %H:%MZ")}  '
                 f'[{i+1}/14]  conf={conf:.0f}')
    ax.legend(fontsize=8, loc='upper right')

    out_path = FRAME_DIR / f'frame_{i+1:02d}_{dt.strftime("%H%M%S")}.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved {out_path.name}')

# ── summary panel: all centroids on one plot ──────────────────────────────
fig, ax = plt.subplots(figsize=(9, 8))

# background: first scan's column-max for spatial reference
col_max_0 = centroids[0][5].max(axis=2)
ax.pcolormesh(x_km, y_km, col_max_0.T,
              cmap=cmap_dbz, norm=norm, shading='auto', alpha=0.4)

# centroid track colored by time
t_vals = np.arange(len(centroids))
cx_all = np.array([c[1] for c in centroids])
cy_all = np.array([c[2] for c in centroids])
cz_all = np.array([c[3] for c in centroids])
cf_all = np.array([c[4] for c in centroids])

ax.plot(cx_all, cy_all, 'w-', lw=1.5, alpha=0.6, zorder=3)
sc = ax.scatter(cx_all, cy_all, c=t_vals, cmap='plasma',
                s=80, zorder=4, edgecolors='white', linewidths=0.5)
plt.colorbar(sc, ax=ax, label='Scan index (0=23:59Z, 13=00:56Z)')

# label every other scan to avoid clutter
for i, (dt, cx, cy, *_) in enumerate(centroids):
    if i % 2 == 0:
        ax.annotate(dt.strftime('%H:%MZ'), (cx, cy),
                    fontsize=7, color='white',
                    xytext=(4, 4), textcoords='offset points')

ax.scatter(-36.0, -2.0, marker='*', c='red', s=200,
           edgecolors='white', zorder=5, label='Bug Fire anchor')
ax.scatter(0, 0, marker='+', c='white', s=150, zorder=5, label='KRGX')

ax.set(aspect='equal',
       xlim=(-40, -25),
       ylim=(-5, 10),
       xlabel='X km (aeqd from KRGX)',
       ylabel='Y km (aeqd from KRGX)',
       title='KRGX Bug Fire — probabilistic centroid track\n'
             '2026-08-09 23:59Z to 2026-08-10 00:56Z')
ax.legend(fontsize=8)
plt.tight_layout()
summary_path = OUT_DIR / 'task17_centroid_summary.png'
plt.savefig(summary_path, dpi=150, bbox_inches='tight')
plt.close()
print(f'\nSaved summary: {summary_path}')
print(f'All frames in: {FRAME_DIR}')
