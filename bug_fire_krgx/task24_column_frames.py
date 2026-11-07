"""
Column tracker overlay on radar — Task 24
Plots column_tracker.json results on per-scan radar backgrounds.
"""
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from pathlib import Path
from datetime import datetime

OUT_DIR    = Path('krgx_hourly')
FRAME_DIR  = OUT_DIR / 'column_frames'
FRAME_DIR.mkdir(exist_ok=True)

ANCHOR_X_KM   = -35.8
ANCHOR_Y_KM   =   1.1
ANCHOR_RAD_KM = 10.0
RADAR_ALT_M   = 2530.0
DBZ_THRESHOLD = 10.0
RHV_MAX       = 0.90

# grid edges — must match krgx_hourly.py exactly
X_MIN, X_MAX = -50_000, 10_000
Y_MIN, Y_MAX = -10_000, 25_000
Z_MAX        =   8_000
DXY, DZ      = 250, 250
x_bins = np.arange(X_MIN, X_MAX + DXY, DXY) / 1e3
y_bins = np.arange(Y_MIN, Y_MAX + DXY, DXY) / 1e3
z_bins = np.arange(0, Z_MAX/1e3 - RADAR_ALT_M/1e3 + DZ/1e3, DZ/1e3)
x_km   = (x_bins[:-1] + x_bins[1:]) / 2
y_km   = (y_bins[:-1] + y_bins[1:]) / 2
z_km_msl = (z_bins[:-1] + z_bins[1:]) / 2 + RADAR_ALT_M/1e3
nx, ny, nz = len(x_km), len(y_km), len(z_km_msl)
XX, YY = np.meshgrid(x_km, y_km, indexing='ij')

# cylinder mask
cyl_mask = ((XX - ANCHOR_X_KM)**2 +
            (YY - ANCHOR_Y_KM)**2) <= ANCHOR_RAD_KM**2
cyl_3d   = cyl_mask[:, :, np.newaxis].repeat(nz, axis=2)

# load column tracker results
with open(OUT_DIR / 'column_tracker.json') as f:
    col_records = json.load(f)
col_by_dt = {r['dt'][:16]: r for r in col_records}  # key by YYYY-MM-DDTHH:MM

# load grid files
grid_files = sorted(OUT_DIR.glob('grids_*.npz'))
norm       = mcolors.Normalize(vmin=0, vmax=35)
cmap_dbz   = plt.cm.jet

for i, gf in enumerate(grid_files):
    stem = gf.stem.replace('grids_', '')
    dt   = datetime.strptime(stem, '%Y%m%d_%H%M%S')

    d   = np.load(gf)
    dbz = d['dbz_smooth'].astype(float)
    rhv = d['rhv_grid'].astype(float)

    # 2D column-max dBZ (full domain)
    col_max = dbz.max(axis=2)

    # smoke mask footprint inside cylinder
    smoke_cyl = (
        cyl_3d &
        np.isfinite(rhv) &
        (rhv < RHV_MAX) &
        (dbz > DBZ_THRESHOLD)
    )

    # burn footprint: lowest 3 z levels with echo inside cylinder
    echo_cyl   = cyl_3d & (dbz > DBZ_THRESHOLD)
    echo_per_z = echo_cyl.sum(axis=(0, 1))
    base_idx   = int(np.argmax(echo_per_z > 0)) if echo_per_z.any() else 0
    burn_slice = echo_cyl[:, :, base_idx:min(base_idx+3, nz)]
    burn_fp    = burn_slice.any(axis=2)   # (nx, ny)

    # smoke footprint column-any inside cylinder
    smoke_fp = smoke_cyl.any(axis=2)

    # get column tracker record for this scan
    dt_key = dt.strftime('%Y-%m-%dT%H:%M')
    rec    = col_by_dt.get(dt_key, {})
    cx     = rec.get('cx')
    cy     = rec.get('cy')
    col_top = rec.get('col_top_msl')
    col_depth = rec.get('col_depth_km')
    burn_area = rec.get('burn_area_km2')

    # ── plot ─────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 8))

    # dBZ background
    pcm = ax.pcolormesh(x_km, y_km, col_max.T,
                        cmap=cmap_dbz, norm=norm, shading='auto')
    plt.colorbar(pcm, ax=ax, label='dBZ (column max)', pad=0.02)

    # smoke mask footprint inside cylinder (purple tint)
    if smoke_fp.any():
        ax.contourf(x_km, y_km, smoke_fp.T.astype(float),
                    levels=[0.5, 1.5], colors=['purple'],
                    alpha=0.20, zorder=2)
        ax.contour(x_km, y_km, smoke_fp.T.astype(float),
                   levels=[0.5], colors=['purple'],
                   linewidths=1.0, linestyles='--', zorder=3)

    # burn footprint (orange tint, lowest gates)
    if burn_fp.any():
        ax.contourf(x_km, y_km, burn_fp.T.astype(float),
                    levels=[0.5, 1.5], colors=['orange'],
                    alpha=0.30, zorder=2)
        ax.contour(x_km, y_km, burn_fp.T.astype(float),
                   levels=[0.5], colors=['darkorange'],
                   linewidths=1.5, linestyles='-', zorder=3)

    # 10 km search cylinder boundary
    theta  = np.linspace(0, 2*np.pi, 360)
    cyl_cx = ANCHOR_X_KM + ANCHOR_RAD_KM * np.cos(theta)
    cyl_cy = ANCHOR_Y_KM + ANCHOR_RAD_KM * np.sin(theta)
    ax.plot(cyl_cx, cyl_cy, 'w-', lw=1.5, alpha=0.7,
            linestyle=':', zorder=4, label='10 km cylinder')

    # GOES anchor
    ax.scatter(ANCHOR_X_KM, ANCHOR_Y_KM,
               marker='*', c='red', s=200,
               edgecolors='white', linewidths=0.5,
               zorder=6, label='GOES anchor')

    # column centroid
    if cx is not None and cy is not None:
        ax.scatter(cx, cy, marker='o', c='yellow', s=180,
                   edgecolors='black', linewidths=1.0,
                   zorder=7, label='Column centroid')
        label_lines = [f'{dt.strftime("%H:%MZ")}']
        if col_top  is not None:
            label_lines.append(f'top={col_top:.2f} km MSL')
        if col_depth is not None:
            label_lines.append(f'depth={col_depth:.2f} km')
        if burn_area is not None:
            label_lines.append(f'burn={burn_area:.1f} km²')
        ax.annotate(
            '\n'.join(label_lines),
            xy=(cx, cy),
            xytext=(10, 10), textcoords='offset points',
            fontsize=8, color='yellow',
            bbox=dict(boxstyle='round,pad=0.3',
                      fc='black', alpha=0.6),
            zorder=8
        )

    # radar site
    ax.scatter(0, 0, marker='+', c='white', s=150,
               linewidths=2, zorder=6, label='KRGX')

    # legend patches for footprints
    smoke_patch = mpatches.Patch(color='purple', alpha=0.4,
                                 label='Smoke mask (RhoHV)')
    burn_patch  = mpatches.Patch(color='orange', alpha=0.5,
                                 label='Burn footprint (low gates)')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles + [smoke_patch, burn_patch],
              fontsize=8, loc='upper right')

    ax.set(
        aspect='equal',
        xlim=(x_km.min(), x_km.max()),
        ylim=(y_km.min(), y_km.max()),
        xlabel='X km (aeqd from KRGX)',
        ylabel='Y km (aeqd from KRGX)',
        title=(f'KRGX  {dt.strftime("%Y-%m-%d %H:%MZ")}  '
               f'[{i+1}/14]\n'
               f'Column: top={col_top:.2f} km  '
               f'depth={col_depth:.2f} km  '
               f'burn={burn_area:.1f} km²'
               if all(v is not None for v in [col_top,col_depth,burn_area])
               else f'KRGX  {dt.strftime("%Y-%m-%d %H:%MZ")}  [{i+1}/14]')
    )

    out_path = FRAME_DIR / f'col_{i+1:02d}_{dt.strftime("%H%M%S")}.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved {out_path.name}  '
          f'top={col_top}  burn={burn_area}')

print(f'\nAll frames saved to {FRAME_DIR}')
