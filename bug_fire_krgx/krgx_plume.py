"""
KRGX plume analysis
- Fetch NEXRAD Level-II from S3
- Grid via histogramdd directly from raw gates
- Gaussian smooth + skeletonize
- 3D plotly: raw gates + gridded isosurface + skeleton + terrain
"""

import io, re, gzip, warnings, time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import boto3
from botocore import UNSIGNED
from botocore.config import Config
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pyart
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
from skimage.morphology import skeletonize
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RegularGridInterpolator
import plotly.graph_objects as go
import requests

warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════
# CONFIG
# ═══════════════════════════════════════════════════════

SITE        = 'KRGX'
RADAR_LAT   = 39.7541
RADAR_LON   = -119.4608
RADAR_ALT_M = 2530.0   # KRGX MSL elevation (m)

T_CENTER = datetime(2026, 8, 10, 0, 8, 4, tzinfo=timezone.utc)
T_BEFORE = timedelta(minutes=0)
T_AFTER  = timedelta(minutes=5)

X_MIN, X_MAX = -50_000, 10_000   # m aeqd
Y_MIN, Y_MAX = -10_000,  25_000
Z_MAX        =   8_000             # m MSL ceiling

DXY, DZ = 250, 250   # m

ELEV_MIN, ELEV_MAX = -0.5, 19.5

DBZ_THRESHOLD = 5
GAUSS_SIGMA   = .7

CACHE_DIR = Path('./krgx_cache')
CACHE_DIR.mkdir(exist_ok=True)
BUCKETS = ['unidata-nexrad-level2', 'noaa-nexrad-level2']

# spectrum width viz
SW_ISO_MIN    =  2.5   # m/s — floor; below this is background turbulence
SW_ISO_MAX    = 12.0   # m/s — ceiling; above this is usually noise/AP
SW_SURF_COUNT =  3
SW_OPACITY    = 0.25

# Anchor derived from GOES-18 FDCC power-weighted fire pixel centroid
# (Task 08b) — independently verified, not derived from radar plume location.
# 'Spot Fire NE' removed: zero independent GOES detection at that location;
# Plume 1 is likely a secondary lobe of the Bug Fire, not a separate source.
FIRE_SOURCES = [
    (39.7631, -119.8788, 'Bug Fire (main)'),
]

# Minimum GOES fire pixels required to call a cluster a confirmed source
GOES_CONFIRMATION_MIN = 1  # set to 0 to disable; currently unused in single-scan

from pyproj import Proj, Transformer

_aeqd  = Proj(proj='aeqd', lat_0=RADAR_LAT, lon_0=RADAR_LON, units='m')
_wgs84 = Proj(proj='latlong', datum='WGS84')
_to_aeqd = Transformer.from_proj(_wgs84, _aeqd, always_xy=True)

fire_source_xy = []   # (x_km, y_km, name) in radar aeqd
for (lat, lon, name) in FIRE_SOURCES:
    x_m, y_m = _to_aeqd.transform(lon, lat)
    fire_source_xy.append((x_m/1e3, y_m/1e3, name))
    print(f"Fire source '{name}': x={x_m/1e3:.1f} km  y={y_m/1e3:.1f} km")

# ═══════════════════════════════════════════════════════
# S3
# ═══════════════════════════════════════════════════════

s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED),
                  region_name='us-east-1')

t_start = T_CENTER - T_BEFORE
t_end   = T_CENTER + T_AFTER

KEY_RE = re.compile(rf'{SITE}(\d{{8}})_(\d{{6}})_V06(\.gz)?$')

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
# TERRAIN  (Open-Elevation API — SRTM 30m)
# ═══════════════════════════════════════════════════════

def fetch_terrain(x_km, y_km, radar_lat, radar_lon):
    from pyproj import Proj
    proj = Proj(proj='aeqd', lat_0=radar_lat, lon_0=radar_lon, datum='WGS84')
    xs = x_km[::4] * 1e3
    ys = y_km[::4] * 1e3
    XX, YY = np.meshgrid(xs, ys)
    lons, lats = proj(XX, YY, inverse=True)
    pts_flat = [{'latitude': float(la), 'longitude': float(lo)}
                for la, lo in zip(lats.ravel(), lons.ravel())]
    elev_flat = []
    for i in range(0, len(pts_flat), 500):
        r = requests.post('https://api.open-elevation.com/api/v1/lookup',
                          json={'locations': pts_flat[i:i+500]}, timeout=30)
        elev_flat.extend([p['elevation'] for p in r.json()['results']])
    elev_coarse = np.array(elev_flat).reshape(XX.shape) / 1e3
    interp = RegularGridInterpolator(
        (y_km[::4], x_km[::4]), elev_coarse, method='linear',
        bounds_error=False, fill_value=None
    )
    YYf, XXf = np.meshgrid(y_km, x_km, indexing='ij')
    return interp((YYf, XXf))   # (ny, nx) km MSL


# ═══════════════════════════════════════════════════════
# FETCH
# ═══════════════════════════════════════════════════════

print('Listing scans...')
scan_list = list_keys(SITE, t_start, t_end)
for _, key, dt in scan_list:
    print(f'  {dt.strftime("%H:%M:%SZ")}  {key.split("/")[-1]}')

bucket, key, fname_dt = scan_list[0]
print(f'\nFetching {key.split("/")[-1]}...')
raw = fetch_raw(bucket, key)

print('Reading radar...')
radar = pyart.io.read_nexrad_archive(
    io.BytesIO(raw), include_fields=[
        'reflectivity', 'spectrum_width',
        'differential_reflectivity', 'cross_correlation_ratio', 'differential_phase'
    ]
)
keep = [
    i for i in range(radar.nsweeps)
    if ELEV_MIN <= float(radar.elevation['data'][radar.get_slice(i)].mean()) <= ELEV_MAX
]
print(f'  {len(keep)} sweeps: '
      f'{[round(float(radar.elevation["data"][radar.get_slice(i)].mean()),1) for i in keep]}')

# ═══════════════════════════════════════════════════════
# RAW GATES  (stored as AGL km)
# ═══════════════════════════════════════════════════════

print('Extracting raw gates...')
gate_x, gate_y, gate_z, gate_v, gate_sw = [], [], [], [], []
gate_zdr, gate_rhv, gate_phidp = [], [], []
for isw in keep:
    sl   = radar.get_slice(isw)
    refl = radar.fields['reflectivity']['data'][sl]
    if np.ma.is_masked(refl): refl = refl.filled(np.nan)
    refl   = refl.astype(float)

    def _extract(field_name, shape):
        f = radar.fields.get(field_name)
        if f is None:
            return np.full(shape, np.nan)
        d = f['data'][sl]
        if np.ma.is_masked(d): d = d.filled(np.nan)
        return d.astype(float)

    sw    = _extract('spectrum_width', refl.shape)
    zdr   = _extract('differential_reflectivity', refl.shape)
    rhv   = _extract('cross_correlation_ratio', refl.shape)
    phidp = _extract('differential_phase', refl.shape)

    gx     = radar.gate_x['data'][sl] / 1e3           # km, aeqd
    gy     = radar.gate_y['data'][sl] / 1e3
    gz_msl = radar.gate_altitude['data'][sl] / 1e3    # km MSL
    gz_agl = gz_msl - RADAR_ALT_M / 1e3               # km AGL
    valid  = (np.isfinite(refl) &
              (gx     >= X_MIN/1e3) & (gx     <= X_MAX/1e3) &
              (gy     >= Y_MIN/1e3) & (gy     <= Y_MAX/1e3) &
              (gz_agl >= 0))                           # above radar site
    gate_x.append(gx[valid].ravel())
    gate_y.append(gy[valid].ravel())
    gate_z.append(gz_agl[valid].ravel())               # AGL
    gate_v.append(refl[valid].ravel())
    gate_sw.append(sw[valid].ravel())
    gate_zdr.append(zdr[valid].ravel())
    gate_rhv.append(rhv[valid].ravel())
    gate_phidp.append(phidp[valid].ravel())

del radar
gw_all     = np.concatenate(gate_sw)
gzdr_all   = np.concatenate(gate_zdr)
grhv_all   = np.concatenate(gate_rhv)
gphidp_all = np.concatenate(gate_phidp)


gx_all = np.concatenate(gate_x)   # km aeqd
gy_all = np.concatenate(gate_y)
gz_all = np.concatenate(gate_z)   # km AGL
gv_all = np.concatenate(gate_v)   # dBZ
print(f'  {len(gx_all)} total gates in domain')

# all points above high threshold first, then fill remainder randomly
HIGH_THRESH = 20
N_PLOT = 50_000

high = (gv_all > HIGH_THRESH)
gx_high, gy_high, gz_high, gv_high = gx_all[high], gy_all[high], gz_all[high], gv_all[high]

remaining = ~high
gx_low, gy_low, gz_low, gv_low = gx_all[remaining], gy_all[remaining], gz_all[remaining], gv_all[remaining]

n_fill = max(0, N_PLOT - len(gx_high))
if len(gx_low) > n_fill:
    idx = np.random.choice(len(gx_low), n_fill, replace=False)
    gx_low, gy_low, gz_low, gv_low = gx_low[idx], gy_low[idx], gz_low[idx], gv_low[idx]

gx_r = np.concatenate([gx_high, gx_low])
gy_r = np.concatenate([gy_high, gy_low])
gz_r = np.concatenate([gz_high, gz_low])
gv_r = np.concatenate([gv_high, gv_low])
gz_r_msl = gz_r + RADAR_ALT_M / 1e3
print(f'  {len(gx_high)} high-dBZ + {len(gx_low)} fill = {len(gx_r)} plot points')

# ═══════════════════════════════════════════════════════
# GRID via histogramdd  (bins in AGL, labels converted to MSL for plotting)
# ═══════════════════════════════════════════════════════

print('Gridding via histogramdd...')
x_bins = np.arange(X_MIN, X_MAX + DXY, DXY) / 1e3
y_bins = np.arange(Y_MIN, Y_MAX + DXY, DXY) / 1e3
z_bins = np.arange(0, Z_MAX/1e3 - RADAR_ALT_M/1e3 + DZ/1e3, DZ/1e3)  # AGL

pts = np.column_stack([gx_all, gy_all, gz_all])



H,      _ = np.histogramdd(pts, bins=[x_bins, y_bins, z_bins], weights=gv_all)
counts, _ = np.histogramdd(pts, bins=[x_bins, y_bins, z_bins])


SW,     _ = np.histogramdd(pts, bins=[x_bins, y_bins, z_bins], weights=gw_all)
sw_grid   = np.where(counts > 0, SW / counts, np.nan)
sw_smooth = gaussian_filter(np.nan_to_num(sw_grid, nan=0.0), sigma=GAUSS_SIGMA)

with np.errstate(invalid='ignore'):
    dbz_grid = np.where(counts > 0, H / counts, np.nan)  # (nx, ny, nz)

# --- polarimetric grids (ZDR, RhoHV, PhiDP) ---
# each field uses its own finite mask + counts, since ZDR/RhoHV/PhiDP can be
# masked independently of reflectivity (avoids biasing the mean low from
# dividing by counts that include gates with no valid polarimetric return)
def _pol_grid(vals_all):
    finite = np.isfinite(vals_all)
    S, _ = np.histogramdd(pts[finite], bins=[x_bins, y_bins, z_bins], weights=vals_all[finite])
    C, _ = np.histogramdd(pts[finite], bins=[x_bins, y_bins, z_bins])
    with np.errstate(invalid='ignore'):
        return np.where(C > 0, S / C, np.nan)

zdr_grid   = _pol_grid(gzdr_all)
rhv_grid   = _pol_grid(grhv_all)
phidp_grid = _pol_grid(gphidp_all)

x_km     = (x_bins[:-1] + x_bins[1:]) / 2    # km aeqd
y_km     = (y_bins[:-1] + y_bins[1:]) / 2
z_km     = (z_bins[:-1] + z_bins[1:]) / 2    # km AGL
z_km_msl = z_km + RADAR_ALT_M / 1e3          # km MSL — used for all plotting

nx, ny, nz = len(x_km), len(y_km), len(z_km)
print(f'  Grid shape: ({nx}, {ny}, {nz})')

ix, iy, iz = np.unravel_index(np.nanargmax(dbz_grid), dbz_grid.shape)
print(f'  Grid max: x={x_km[ix]:.1f} y={y_km[iy]:.1f} '
      f'z={z_km_msl[iz]:.2f} km MSL  val={np.nanmax(dbz_grid):.1f} dBZ')

# ═══════════════════════════════════════════════════════
# GAUSSIAN SMOOTH + SKELETONIZE
# ═══════════════════════════════════════════════════════

print('Skeletonizing...')
dbz_smooth = gaussian_filter(np.nan_to_num(dbz_grid, nan=0.0), sigma=GAUSS_SIGMA)
binary_3d  = dbz_smooth > 15

# ═══════════════════════════════════════════════════════
# POLARIMETRIC SMOKE MASK
# ═══════════════════════════════════════════════════════
# Requires: valid dBZ, low RhoHV (irregular particles), ZDR in smoke range
RHV_MAX = 0.90   # low RhoHV = irregular/tumbling particles (smoke/ash)

smoke_mask_3d = (
    (dbz_smooth > DBZ_THRESHOLD) &       # non-trivial echo
    np.isfinite(rhv_grid) &              # RhoHV actually observed here
    (rhv_grid < RHV_MAX)
)
print(f"Smoke mask: {smoke_mask_3d.sum()} voxels "
      f"({100*smoke_mask_3d.mean():.2f}% of grid)")
skel_3d    = skeletonize(binary_3d)
sv_x, sv_y, sv_z = np.where(skel_3d)
skel_xkm = x_km[sv_x]
skel_ykm = y_km[sv_y]
skel_zkm = z_km_msl[sv_z]   # MSL for plotting
print(f'  {skel_3d.sum()} skeleton voxels')

# ═══════════════════════════════════════════════════════
# PCA per cluster
# ═══════════════════════════════════════════════════════

xx3, yy3, zz3 = np.meshgrid(x_km, y_km, z_km_msl, indexing='ij')  # all MSL
pts_bin = np.column_stack([xx3[binary_3d], yy3[binary_3d], zz3[binary_3d]])

plume_results = []   # one dict per DBSCAN cluster; reused by the hourly pipeline (Task 09)

if len(pts_bin) > 10:
    db      = DBSCAN(eps=3, min_samples=5).fit(pts_bin)
    labels  = db.labels_
    n_clust = len(set(labels)) - (1 if -1 in labels else 0)
    print(f'\n{n_clust} plume clusters')
    for cid in range(n_clust):
        cpts = pts_bin[labels == cid]
        pca  = PCA(n_components=3).fit(cpts)
        c    = cpts.mean(axis=0)
        d    = pca.components_[0]
        elev = np.degrees(np.arctan2(d[2], np.sqrt(d[0]**2 + d[1]**2)))
        hor  = np.degrees(np.arctan2(d[1], d[0]))
        cx, cy = c[0], c[1]

        # Nearest fire source
        dists = [np.sqrt((cx - sx)**2 + (cy - sy)**2)
                 for (sx, sy, _) in fire_source_xy]
        nearest_idx  = int(np.argmin(dists))
        nearest_name = fire_source_xy[nearest_idx][2]
        nearest_dist = dists[nearest_idx]

        # A cluster is a confirmed source lobe if it has a named anchor
        # and dist is plausible (< 20 km — smoke advects, so generous bound).
        # Clusters with no nearby anchor or dist > 20 km → unconfirmed.
        DIST_THRESHOLD_KM = 20.0
        if nearest_dist <= DIST_THRESHOLD_KM:
            source_label     = nearest_name
            source_confirmed = True
        else:
            source_label     = f'Unconfirmed (nearest: {nearest_name} {nearest_dist:.1f} km)'
            source_confirmed = False

        print(f'  Plume {cid}: {len(cpts)} pts  '
              f'x={c[0]:.1f} y={c[1]:.1f} z={c[2]:.2f} km MSL  '
              f'elev={elev:.1f}°  horiz={hor:.1f}°')
        print(f"  Plume {cid} -> '{source_label}' "
              f"({'confirmed' if source_confirmed else 'UNCONFIRMED'}, "
              f"dist={nearest_dist:.1f} km)")

        plume_results.append({
            'cluster_id': cid,
            'n_pts': len(cpts),
            'x_km': cx, 'y_km': cy, 'z_km_msl': c[2],
            'elev_deg': elev, 'horiz_deg': hor,
            'source_idx': nearest_idx,
            'source_name': nearest_name,
            'source_confirmed': source_confirmed,
            'source_dist_km': round(nearest_dist, 2),
        })

    # per-plume overlap with the polarimetric smoke mask
    # (labels are per-voxel over pts_bin = xx3[binary_3d] etc, in the same
    # boolean-index order as binary_3d, so we can scatter them back into a
    # full-grid label array to compare against smoke_mask_3d voxel-for-voxel)
    label_grid = np.full(binary_3d.shape, -1, dtype=int)
    label_grid[binary_3d] = labels
    for cid in range(n_clust):
        plume_mask = (label_grid == cid)
        plume_vox  = plume_mask.sum()
        overlap    = (smoke_mask_3d & plume_mask).sum()
        print(f'  Plume {cid}: {plume_vox} vox total, '
              f'{overlap} ({100*overlap/max(plume_vox,1):.0f}%) in smoke mask')

# ═══════════════════════════════════════════════════════
# TERRAIN
# ═══════════════════════════════════════════════════════

terrain_cache = CACHE_DIR / 'terrain.npy'
if terrain_cache.exists():
    print('Loading cached terrain...')
    terrain_km = np.load(terrain_cache)
else:
    print('Fetching terrain (SRTM via Open-Elevation)...')
    terrain_km = fetch_terrain(x_km, y_km, RADAR_LAT, RADAR_LON)
    np.save(terrain_cache, terrain_km)
    print(f'terrain shape: {terrain_km.shape}  expected ({len(y_km)}, {len(x_km)})')

YY_t, XX_t = np.meshgrid(y_km, x_km)   # (nx, ny) — x varies along rows

# ═══════════════════════════════════════════════════════
# 3D PLOTLY
# ═══════════════════════════════════════════════════════

print('Building 3D plot...')
vf_grid = np.where(np.isfinite(dbz_smooth), dbz_smooth, np.nan)

fig3d = go.Figure()

# # terrain
# fig3d.add_trace(go.Surface(
#     x=XX_t, y=YY_t, z=terrain_km.T,   # transpose terrain to (nx, ny)
#     colorscale='earth',
#     showscale=False,
#     opacity=0.7,
#     name='terrain',
# ))

# raw gates (MSL) — only above threshold
thresh_mask = gv_r > DBZ_THRESHOLD
fig3d.add_trace(go.Scatter3d(
    x=gx_r[thresh_mask], y=gy_r[thresh_mask], z=gz_r_msl[thresh_mask],
    mode='markers',
    name='raw gates',
    marker=dict(
        symbol='circle', size=5,
        color=gv_r[thresh_mask], colorscale='Jet',
        cmin=0, cmax=30, opacity=0.3,
        colorbar=dict(title='dBZ', x=1.05),
    )
))

# gridded isosurface (MSL via zz3)
fig3d.add_trace(go.Isosurface(
    x=xx3.ravel(), y=yy3.ravel(), z=zz3.ravel(),
    value=vf_grid.ravel(),
    isomin=DBZ_THRESHOLD,
    isomax=30,
    surface_count=10,
    colorscale='jet',
    opacity=0.4,
    showscale=False,
    caps=dict(x_show=False, y_show=False, z_show=False),
    name='gridded',
))

# skeleton (MSL)
skel_dbz = dbz_smooth[sv_x, sv_y, sv_z]
fig3d.add_trace(go.Scatter3d(
    x=skel_xkm, y=skel_ykm, z=skel_zkm,
    mode='markers',
    name='skeleton',
    marker=dict(
        symbol='diamond', size=5,
        color=skel_dbz, colorscale='jet',
        opacity=1.0,
    )
))

# radar site marker
fig3d.add_trace(go.Scatter3d(
    x=[0], y=[0], z=[RADAR_ALT_M/1e3],
    mode='markers+text',
    name='KRGX',
    text=['KRGX'],
    textposition='top center',
    marker=dict(symbol='diamond', size=8, color='white'),
))

z_plot_min = 1.0   # km MSL — below terrain   
#z_plot_min = RADAR_ALT_M / 1e3
z_plot_max = Z_MAX / 1e3

fig3d.update_layout(
    scene=dict(
        xaxis=dict(title='X (km E+)'),
        yaxis=dict(title='Y (km N+)'),
        zaxis=dict(title='Z MSL (km)', range=[z_plot_min, z_plot_max]),
        aspectmode='manual',
        aspectratio=dict(
            x=(X_MAX - X_MIN) / (Y_MAX - Y_MIN),
            y=1.0,
            z=(z_plot_max - z_plot_min) / ((Y_MAX - Y_MIN)/1e3) * 3,
        ),
        camera=dict(eye=dict(x=1.4, y=-1.4, z=0.9)),
    ),
    width=1100, height=800,
    title=f'KRGX  {fname_dt.strftime("%H:%M:%S Z")}  |  KRGX {RADAR_ALT_M/1e3:.2f} km MSL',
    legend=dict(x=0, y=1),
)

vf_sw = np.where(np.isfinite(sw_smooth), sw_smooth, np.nan)

# RhoHV smoke-mask volume — show only where smoke_mask_3d is True
rhv_plot = np.where(smoke_mask_3d, rhv_grid, np.nan)
fig3d.add_trace(go.Volume(
    x=xx3.ravel(), y=yy3.ravel(), z=zz3.ravel(),
    value=rhv_plot.ravel(),
    isomin=0.50,
    isomax=0.90,
    surface_count=5,
    colorscale='Purples_r',   # dark purple = low RhoHV = most smoke-like
    opacity=0.15,
    caps=dict(x_show=False, y_show=False, z_show=False),
    name='Smoke mask (ρHV)',
    showscale=True,
    colorbar=dict(title='ρHV', x=1.2, len=0.5),
))

fig3d.add_trace(go.Isosurface(
    x=xx3.ravel(), y=yy3.ravel(), z=zz3.ravel(),
    value=vf_sw.ravel(),
    isomin=SW_ISO_MIN,
    isomax=SW_ISO_MAX,
    surface_count=SW_SURF_COUNT,
    colorscale=[
        [0.0, 'yellow'],
        [0.5, 'orange'],
        [1.0, 'orangered'],
    ],
    opacity=SW_OPACITY,
    flatshading=True,
    lighting=dict(
        ambient=0.4,
        diffuse=0.8,
        specular=0.1,
        roughness=0.9,
        fresnel=0.0,
    ),
    lightposition=dict(x=100, y=100, z=1000),

    surface=dict(
        show=True,
        count=SW_SURF_COUNT,
    ),
    caps=dict(x_show=False, y_show=False, z_show=False),
    contour=dict(
        show=True,
        color='black',
        width=1,
    ),
    showscale=True,
    colorbar=dict(
        title='SW (m/s)',
        x=1.1, thickness=16, len=0.65,
        tickvals=np.arange(SW_ISO_MIN, SW_ISO_MAX+1, 2).tolist(),
    ),
    name='spectrum width',
))




fig3d.add_trace(go.Surface(
    x=XX_t, y=YY_t, z=terrain_km,
    colorscale='earth',
    showscale=False,
    opacity=0.7,
    name='terrain',
))

fig3d.show()
fig3d.write_html('krgx_3d_plume.html', include_plotlyjs='cdn')
print('Saved krgx_3d_plume.html')
