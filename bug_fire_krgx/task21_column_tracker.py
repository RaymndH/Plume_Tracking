"""
Rising column tracker — Task 21
Extracts vertical structure of the fire injection column within a
search cylinder around the Bug Fire GOES anchor across 14 scans.
"""
import io, json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
import pyart, boto3
from botocore import UNSIGNED
from botocore.config import Config

OUT_DIR   = Path('krgx_hourly')
CACHE_DIR = Path('krgx_cache')

# ── fire anchor (GOES power-weighted centroid, Task 08b) ──────────────────
ANCHOR_X_KM  = -35.8
ANCHOR_Y_KM  =   1.1
ANCHOR_RAD_KM = 10.0   # search cylinder radius

# ── radar config (must match krgx_hourly.py) ─────────────────────────────
RADAR_LAT   = 39.7541
RADAR_LON   = -119.4608
RADAR_ALT_M = 2530.0
X_MIN, X_MAX = -50_000, 10_000
Y_MIN, Y_MAX = -10_000, 25_000
Z_MAX        =   8_000
DXY, DZ      = 250, 250
ELEV_MIN, ELEV_MAX = -0.5, 19.5
DBZ_THRESHOLD = 10.0
RHV_MAX       = 0.90

# ── grid edges (must match krgx_hourly.py exactly) ───────────────────────
x_bins = np.arange(X_MIN, X_MAX + DXY, DXY) / 1e3
y_bins = np.arange(Y_MIN, Y_MAX + DXY, DXY) / 1e3
z_bins = np.arange(0, Z_MAX/1e3 - RADAR_ALT_M/1e3 + DZ/1e3, DZ/1e3)
x_km   = (x_bins[:-1] + x_bins[1:]) / 2
y_km   = (y_bins[:-1] + y_bins[1:]) / 2
z_km   = (z_bins[:-1] + z_bins[1:]) / 2
z_km_msl = z_km + RADAR_ALT_M / 1e3

nx, ny, nz = len(x_km), len(y_km), len(z_km)
XX, YY = np.meshgrid(x_km, y_km, indexing='ij')   # (nx, ny)

# ── cylinder mask (2D, applied at all z levels) ───────────────────────────
dist2_xy  = (XX - ANCHOR_X_KM)**2 + (YY - ANCHOR_Y_KM)**2
cyl_mask  = dist2_xy <= ANCHOR_RAD_KM**2   # (nx, ny) bool

# expand to 3D
cyl_3d = cyl_mask[:, :, np.newaxis].repeat(nz, axis=2)   # (nx, ny, nz)

# ── load Vr per scan (need raw Level-II for velocity) ────────────────────
def load_vr_grid(cache_path):
    """Load velocity field and grid it the same way as dbz/rhv."""
    raw   = cache_path.read_bytes()
    radar = pyart.io.read_nexrad_archive(
        io.BytesIO(raw),
        include_fields=['reflectivity', 'velocity', 'cross_correlation_ratio']
    )
    keep = [
        i for i in range(radar.nsweeps)
        if ELEV_MIN <= float(
            radar.elevation['data'][radar.get_slice(i)].mean()
        ) <= ELEV_MAX
    ]
    gate_x, gate_y, gate_z = [], [], []
    gate_vr, gate_rhv = [], []

    for isw in keep:
        sl     = radar.get_slice(isw)
        refl   = radar.fields['reflectivity']['data'][sl]
        vr_f   = radar.fields.get('velocity')
        rhv_f  = radar.fields.get('cross_correlation_ratio')

        gx     = radar.gate_x['data'][sl] / 1e3
        gy     = radar.gate_y['data'][sl] / 1e3
        gz_agl = radar.gate_altitude['data'][sl]/1e3 - RADAR_ALT_M/1e3

        valid = (
            ~np.ma.getmaskarray(refl) &
            (gx >= X_MIN/1e3) & (gx <= X_MAX/1e3) &
            (gy >= Y_MIN/1e3) & (gy <= Y_MAX/1e3) &
            (gz_agl >= 0)
        )
        gate_x.append(gx[valid].ravel())
        gate_y.append(gy[valid].ravel())
        gate_z.append(gz_agl[valid].ravel())

        def _extract(field):
            if field is None:
                return np.full(valid.sum(), np.nan)
            d = field['data'][sl]
            return np.where(np.ma.getmaskarray(d), np.nan,
                            d.astype(float))[valid].ravel()

        gate_vr.append(_extract(vr_f))
        gate_rhv.append(_extract(rhv_f))

    del radar
    if not gate_x:
        return None, None

    pts    = np.column_stack([
        np.concatenate(gate_x),
        np.concatenate(gate_y),
        np.concatenate(gate_z)
    ])
    gvr    = np.concatenate(gate_vr)
    grhv   = np.concatenate(gate_rhv)

    def _pol_grid(vals):
        finite = np.isfinite(vals)
        if finite.sum() == 0:
            return np.full((nx,ny,nz), np.nan)
        s, _ = np.histogramdd(pts[finite], bins=[x_bins,y_bins,z_bins],
                              weights=vals[finite])
        c, _ = np.histogramdd(pts[finite], bins=[x_bins,y_bins,z_bins])
        return np.where(c > 0, s/c, np.nan)

    return _pol_grid(gvr), _pol_grid(grhv)

# ── main loop ─────────────────────────────────────────────────────────────
grid_files  = sorted(OUT_DIR.glob('grids_*.npz'))
cache_files = sorted(CACHE_DIR.glob('KRGX20260810_*'))
# also include the 23:59 scan from 20260809
cache_files += sorted(CACHE_DIR.glob('KRGX20260809_*'))
cache_files  = sorted(set(cache_files))

# build dt→cache_path lookup
cache_lookup = {}
for cf in cache_files:
    stem = cf.stem   # e.g. KRGX20260810_000404_V06
    parts = stem.split('_')
    if len(parts) >= 2:
        try:
            dt = datetime.strptime(parts[0][4:] + parts[1], '%Y%m%d%H%M%S')
            cache_lookup[dt] = cf
        except: pass

records = []

prev_cx, prev_cy = None, None   # seeded from first scan's raw centroid
MAX_DISP_KM = 2.0   # max allowed centroid displacement per scan

for gf in grid_files:
    stem = gf.stem.replace('grids_', '')
    dt   = datetime.strptime(stem, '%Y%m%d_%H%M%S')

    d         = np.load(gf)
    dbz       = d['dbz_smooth'].astype(float)
    rhv       = d['rhv_grid'].astype(float)

    # ── restrict to cylinder ─────────────────────────────────────────────
    dbz_cyl = np.where(cyl_3d, dbz, np.nan)
    rhv_cyl = np.where(cyl_3d, rhv, np.nan)

    # smoke mask inside cylinder
    smoke_cyl = (
        np.isfinite(rhv_cyl) &
        (rhv_cyl < RHV_MAX) &
        (dbz_cyl > DBZ_THRESHOLD)
    )
    echo_cyl = dbz_cyl > DBZ_THRESHOLD

    # ── vertical profile inside cylinder ─────────────────────────────────
    # count echo voxels per z level
    echo_per_z  = echo_cyl.sum(axis=(0,1))    # (nz,)
    smoke_per_z = smoke_cyl.sum(axis=(0,1))   # (nz,)

    # column base: lowest z with any echo
    base_idx = np.argmax(echo_per_z > 0) if echo_per_z.any() else None
    # column top: highest z with any echo
    top_idx  = (len(echo_per_z) - 1 - np.argmax(
        echo_per_z[::-1] > 0)) if echo_per_z.any() else None

    col_base_msl = float(z_km_msl[base_idx]) if base_idx is not None else np.nan
    col_top_msl  = float(z_km_msl[top_idx])  if top_idx  is not None else np.nan
    col_depth_km = (col_top_msl - col_base_msl) if np.isfinite(col_base_msl) else np.nan

    # column centroid z (echo-weighted)
    dbz_col_max = dbz_cyl.max(axis=2)   # (nx,ny) column max inside cylinder
    total_echo  = echo_per_z.sum()
    if total_echo > 0:
        # x/y centroid from column-max footprint (raw dBZ — used as fallback)
        w2d = np.where(np.isfinite(dbz_col_max) & (dbz_col_max > DBZ_THRESHOLD),
                       dbz_col_max, 0.0)
        # Weight by smoke-masked dBZ only — low-RhoHV gates near injection
        # point should dominate over advected plume mass
        smoke_col_max = np.where(
            smoke_cyl.any(axis=2),
            np.where(np.isfinite(dbz_col_max), dbz_col_max, 0.0),
            0.0
        )
        if smoke_col_max.sum() > 0:
            cx = float(np.average(x_km, weights=smoke_col_max.sum(axis=1)))
            cy = float(np.average(y_km, weights=smoke_col_max.sum(axis=0)))
        elif w2d.sum() > 0:
            # fallback to raw dBZ if no smoke mask coverage
            cx = float(np.average(x_km, weights=w2d.sum(axis=1)))
            cy = float(np.average(y_km, weights=w2d.sum(axis=0)))
        else:
            cx = cy = np.nan

        smoke_per_z_w = np.where(smoke_per_z > 0, smoke_per_z,
                                 echo_per_z * 0.1)   # downweight non-smoke
        cz = float(np.average(z_km_msl, weights=smoke_per_z_w))
    else:
        cx = cy = cz = np.nan

    # ── scan-to-scan displacement cap (stabilize centroid track) ─────────
    if np.isfinite(cx) and np.isfinite(cy):
        if prev_cx is None:
            # First scan — accept raw centroid unconditionally as seed
            prev_cx, prev_cy = cx, cy
            print(f"    [SEED] {dt.strftime('%H:%MZ')}: "
                  f"seeded at ({cx:.1f}, {cy:.1f})")
        else:
            disp = np.sqrt((cx - prev_cx)**2 + (cy - prev_cy)**2)
            if disp > MAX_DISP_KM:
                print(f"    [CAP] {dt.strftime('%H:%MZ')}: "
                      f"jump {disp:.1f} km > cap — holding "
                      f"({prev_cx:.1f}, {prev_cy:.1f})")
                cx, cy = prev_cx, prev_cy
            else:
                prev_cx, prev_cy = cx, cy
    else:
        if prev_cx is not None:
            cx, cy = prev_cx, prev_cy   # hold last known

    # ── burning region footprint (lowest 3 z levels with echo) ───────────
    # horizontal extent of echo near ground level
    if base_idx is not None:
        low_z_slice = echo_cyl[:, :, base_idx:min(base_idx+3, nz)]
        burn_footprint = low_z_slice.any(axis=2)   # (nx,ny)
        burn_area_km2  = float(burn_footprint.sum() * (DXY/1e3)**2)
        burn_cx = float(np.average(x_km, weights=burn_footprint.sum(axis=1))
                        ) if burn_footprint.any() else np.nan
        burn_cy = float(np.average(y_km, weights=burn_footprint.sum(axis=0))
                        ) if burn_footprint.any() else np.nan
    else:
        burn_footprint = None
        burn_area_km2  = np.nan
        burn_cx = burn_cy = np.nan

    # ── Vr / falling column ───────────────────────────────────────────────
    # Find nearest cache file
    best_cf = min(cache_lookup.items(),
                  key=lambda kv: abs((kv[0]-dt).total_seconds()),
                  default=(None, None))
    falling_flag  = False
    vr_top_ms     = np.nan
    vr_mean_cyl   = np.nan

    if best_cf[1] is not None and abs(
            (best_cf[0]-dt).total_seconds()) < 300:
        vr_grid, _ = load_vr_grid(best_cf[1])
        if vr_grid is not None:
            vr_cyl = np.where(cyl_3d, vr_grid, np.nan)
            # Mean Vr in upper half of column (positive = away from radar)
            if top_idx is not None and base_idx is not None:
                mid_idx  = (base_idx + top_idx) // 2
                vr_upper = vr_cyl[:, :, mid_idx:top_idx+1]
                vr_top_ms = float(np.nanmean(vr_upper))
            vr_mean_cyl = float(np.nanmean(vr_cyl))
            # Falling column flag: Vr at column top is strongly negative
            # (toward radar = downward motion at high elevation angles)
            # Use conservative threshold: mean upper Vr < -3 m/s
            if np.isfinite(vr_top_ms) and vr_top_ms < -3.0:
                falling_flag = True

    rec = dict(
        dt=dt.isoformat(),
        cx=round(cx, 2) if np.isfinite(cx) else None,
        cy=round(cy, 2) if np.isfinite(cy) else None,
        cz=round(cz, 2) if np.isfinite(cz) else None,
        col_top_msl=round(col_top_msl, 2)   if np.isfinite(col_top_msl)  else None,
        col_depth_km=round(col_depth_km, 2) if np.isfinite(col_depth_km) else None,
        burn_area_km2=round(burn_area_km2, 2) if np.isfinite(burn_area_km2) else None,
        burn_cx=round(burn_cx, 2) if np.isfinite(burn_cx) else None,
        burn_cy=round(burn_cy, 2) if np.isfinite(burn_cy) else None,
        vr_top_ms=round(vr_top_ms, 2) if np.isfinite(vr_top_ms) else None,
        vr_mean_cyl=round(vr_mean_cyl, 2) if np.isfinite(vr_mean_cyl) else None,
        falling_flag=falling_flag,
        echo_vox_total=int(total_echo),
    )
    records.append(rec)

    print(f"  {dt.strftime('%H:%MZ')}  "
          f"top={col_top_msl:.2f}  "
          f"depth={col_depth_km:.2f} km  "
          f"burn={burn_area_km2:.1f} km²")

# ── save results ──────────────────────────────────────────────────────────
out_json = OUT_DIR / 'column_tracker.json'
with open(out_json, 'w') as f:
    json.dump(records, f, indent=2)
print(f'\nSaved {out_json}')

# ── plots ─────────────────────────────────────────────────────────────────
t_min = np.array([(datetime.fromisoformat(r['dt']) -
                   datetime.fromisoformat(records[0]['dt'])
                   ).total_seconds()/60 for r in records])

def _arr(key):
    return np.array([r[key] if r[key] is not None else np.nan
                     for r in records])

fig, axes = plt.subplots(2, 2, figsize=(12, 9))

# 1. Column top + centroid z
ax = axes[0, 0]
ax.plot(t_min, _arr('col_top_msl'), 's-', color='firebrick', label='Top')
ax.plot(t_min, _arr('cz'),          '^-', color='darkorange', label='Centroid z')
ax.axhline(RADAR_ALT_M/1e3, color='gray', lw=1, linestyle=':',
           label=f'Radar floor ({RADAR_ALT_M/1e3:.2f} km MSL)')
ax.set(xlabel='Minutes', ylabel='km MSL',
       title='Column top + centroid height\n(base not shown — beam floor artifact)')
ax.legend(fontsize=8)

# 2. Column depth
ax = axes[0, 1]
ax.plot(t_min, _arr('col_depth_km'), 'o-', color='purple')
ax.set(xlabel='Minutes', ylabel='km', title='Column depth (top - base)')

# 3. Burning region area
ax = axes[1, 0]
ax.plot(t_min, _arr('burn_area_km2'), 'o-', color='darkorange')
ax.set(xlabel='Minutes', ylabel='km²', title='Radar-estimated burning area')

# 4. Column centroid x/y drift (stabilized, with 2 km scale bar)
ax = axes[1, 1]
sc = ax.scatter(_arr('cx'), _arr('cy'),
                c=t_min, cmap='plasma', s=80, zorder=4,
                edgecolors='white', linewidths=0.5)
ax.plot(_arr('cx'), _arr('cy'), 'k-', lw=1.0, zorder=3)
for i, lbl in enumerate(
        [datetime.fromisoformat(r['dt']).strftime('%H:%M')
         for r in records]):
    if i % 3 == 0:   # label every 3rd scan
        ax.annotate(lbl, (_arr('cx')[i], _arr('cy')[i]),
                    fontsize=7, xytext=(4,4),
                    textcoords='offset points')
# 2 km scale bar
sb_x = ANCHOR_X_KM + 2.0
sb_y = ANCHOR_Y_KM - 3.0
ax.annotate('', xy=(sb_x + 2, sb_y), xytext=(sb_x, sb_y),
            arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
ax.text(sb_x + 1, sb_y - 0.3, '2 km', ha='center', fontsize=8)
ax.scatter(ANCHOR_X_KM, ANCHOR_Y_KM, marker='*', c='red',
           s=200, zorder=5, label='GOES anchor')
plt.colorbar(sc, ax=ax, label='Minutes from first scan')
ax.set(xlabel='X km', ylabel='Y km',
       title='Column centroid drift (2 km cap applied)',
       aspect='equal')
ax.legend(fontsize=8)

plt.suptitle('KRGX Bug Fire — Column tracker\n'
             '2026-08-09 23:59Z to 2026-08-10 00:56Z', fontsize=12)
plt.tight_layout()
out_png = OUT_DIR / 'task21_column_tracker.png'
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print(f'Saved {out_png}')
