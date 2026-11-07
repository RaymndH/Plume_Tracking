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

SITE        = 'KRGX'
RADAR_LAT   = 39.7541
RADAR_LON   = -119.4608
RADAR_ALT_M = 2530.0   # KRGX MSL elevation (m)

# ── CONFIG ──────────────────────────────────────────────
T_CENTER = datetime(2026, 8, 9, 23, 59, 57, tzinfo=timezone.utc)
T_BEFORE = timedelta(minutes=0)
T_AFTER  = timedelta(minutes=60)



X_MIN, X_MAX = -50_000,  10_000   # m aeqd
Y_MIN, Y_MAX = -10_000,  25_000
Z_MAX        =   8_000             # m MSL ceiling

DXY, DZ = 250, 250   # m

ELEV_MIN, ELEV_MAX = -0.5, 19.5

DBZ_THRESHOLD = 5
GAUSS_SIGMA   = .7

CACHE_DIR = Path('./krgx_cache')
CACHE_DIR.mkdir(exist_ok=True)
BUCKETS = ['unidata-nexrad-level2', 'noaa-nexrad-level2']

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


# ── LIST ALL SCANS ──────────────────────────────────────
print('Listing scans...')
scan_list = list_keys(SITE, t_start, t_end)
print(f'{len(scan_list)} scans found')

# ── PROCESS EACH SCAN ───────────────────────────────────
results = []   # list of dicts, one per scan

for si, (bucket, key, fname_dt) in enumerate(scan_list):
    print(f'\n[{si+1}/{len(scan_list)}] {fname_dt.strftime("%H:%M:%SZ")}')

    # fetch
    try:
        raw = fetch_raw(bucket, key)
    except Exception as e:
        print(f'  fetch error: {e}')
        continue

    # read + extract gates
    radar = pyart.io.read_nexrad_archive(
        io.BytesIO(raw), include_fields=['reflectivity']
    )
    keep = [
        i for i in range(radar.nsweeps)
        if ELEV_MIN <= float(radar.elevation['data'][radar.get_slice(i)].mean()) <= ELEV_MAX
    ]

    gate_x, gate_y, gate_z, gate_v = [], [], [], []
    for sw in keep:
        sl   = radar.get_slice(sw)
        refl = radar.fields['reflectivity']['data'][sl]
        if np.ma.is_masked(refl): refl = refl.filled(np.nan)
        refl   = refl.astype(float)
        gx     = radar.gate_x['data'][sl] / 1e3
        gy     = radar.gate_y['data'][sl] / 1e3
        gz_agl = radar.gate_altitude['data'][sl] / 1e3 - RADAR_ALT_M / 1e3
        valid  = (np.isfinite(refl) &
                  (gx >= X_MIN/1e3) & (gx <= X_MAX/1e3) &
                  (gy >= Y_MIN/1e3) & (gy <= Y_MAX/1e3) &
                  (gz_agl >= 0))
        gate_x.append(gx[valid].ravel())
        gate_y.append(gy[valid].ravel())
        gate_z.append(gz_agl[valid].ravel())
        gate_v.append(refl[valid].ravel())
    del radar

    if not gate_x:
        print('  no gates in domain')
        continue

    gx_all = np.concatenate(gate_x)
    gy_all = np.concatenate(gate_y)
    gz_all = np.concatenate(gate_z)
    gv_all = np.concatenate(gate_v)

    # grid via histogramdd
    x_bins = np.arange(X_MIN, X_MAX + DXY, DXY) / 1e3
    y_bins = np.arange(Y_MIN, Y_MAX + DXY, DXY) / 1e3
    z_bins = np.arange(0, Z_MAX/1e3 - RADAR_ALT_M/1e3 + DZ/1e3, DZ/1e3)

    pts = np.column_stack([gx_all, gy_all, gz_all])
    H,      _ = np.histogramdd(pts, bins=[x_bins, y_bins, z_bins], weights=gv_all)
    counts, _ = np.histogramdd(pts, bins=[x_bins, y_bins, z_bins])
    with np.errstate(invalid='ignore'):
        dbz_grid = np.where(counts > 0, H / counts, np.nan)

    x_km     = (x_bins[:-1] + x_bins[1:]) / 2
    y_km     = (y_bins[:-1] + y_bins[1:]) / 2
    z_km     = (z_bins[:-1] + z_bins[1:]) / 2
    z_km_msl = z_km + RADAR_ALT_M / 1e3

    # smooth + binary
    dbz_smooth = gaussian_filter(np.nan_to_num(dbz_grid, nan=0.0), sigma=GAUSS_SIGMA)
    binary_3d  = dbz_smooth > DBZ_THRESHOLD

    if binary_3d.sum() < 10:
        print('  insufficient echo coverage')
        continue

    # peak detection
    col_max_2d  = dbz_smooth.max(axis=2)
    global_max  = dbz_smooth.max()
    MIN_PEAK_VAL = global_max * 0.65
    MIN_NEIGHBOR_VOXELS = round((2000/DXY)**2 * 0.3)

    coords_2d = peak_local_max(col_max_2d, min_distance=5, threshold_abs=8)

    real_peaks = []
    for px, py in coords_2d:
        peak_val = col_max_2d[px, py]
        r_vox = round(2000 / DXY)
        x0, x1 = max(0, px-r_vox), min(col_max_2d.shape[0], px+r_vox)
        y0, y1 = max(0, py-r_vox), min(col_max_2d.shape[1], py+r_vox)
        n_above = (col_max_2d[x0:x1, y0:y1] > DBZ_THRESHOLD).sum()
        if peak_val >= MIN_PEAK_VAL and n_above >= MIN_NEIGHBOR_VOXELS:
            real_peaks.sort(key=lambda x: -x[2])
            real_peaks.append((px, py, peak_val))

    real_peaks.sort(key=lambda x: -x[2])
    chosen = [(px, py) for px, py, _ in real_peaks]

    if not chosen:
        print('  no real peaks found')
        continue

    # watershed
    markers = np.zeros_like(dbz_smooth, dtype=int)
    for i, (px, py) in enumerate(chosen):
        pz = np.argmax(dbz_smooth[px, py, :])
        markers[px, py, pz] = i + 1

    labels_ws = watershed(-dbz_smooth, markers, mask=binary_3d)

    # saddle-based merging
    from itertools import combinations
    from scipy.ndimage import binary_dilation

    merge_map = {i+1: i+1 for i in range(len(chosen))}  # rid → group

    for r1, r2 in combinations(range(1, len(chosen)+1), 2):
        reg1 = labels_ws == r1
        reg2 = labels_ws == r2
        border = binary_dilation(reg1) & binary_dilation(reg2) & binary_3d
        if border.sum() == 0: continue
        saddle   = dbz_smooth[border].max()
        peak1    = dbz_smooth[reg1].max()
        peak2    = dbz_smooth[reg2].max()
        drop_min = min(peak1 - saddle, peak2 - saddle)
        if drop_min <= 3.0:
            # merge r2 into r1's group
            old = merge_map[r2]
            new = merge_map[r1]
            for k in merge_map:
                if merge_map[k] == old:
                    merge_map[k] = new

    # collect merged regions
    groups = {}
    for rid, gid in merge_map.items():
        region = labels_ws == rid
        groups.setdefault(gid, np.zeros_like(binary_3d, dtype=bool))
        groups[gid] |= region

    scan_result = {'dt': fname_dt, 'hotspots': []}
    for gid, region in groups.items():
        if region.sum() == 0: continue
        ix, iy, iz = np.where(region)
        hs = {
            'x_mean': float(x_km[ix].mean()),
            'y_mean': float(y_km[iy].mean()),
            'z_mean': float(z_km_msl[iz].mean()),
            'z_max':  float(z_km_msl[iz].max()),
            'dbz_max': float(dbz_smooth[region].max()),
            'voxels': int(region.sum()),
        }
        scan_result['hotspots'].append(hs)
        print(f'  Hotspot: x={hs["x_mean"]:.1f} y={hs["y_mean"]:.1f} '
              f'z={hs["z_mean"]:.2f} km MSL  '
              f'top={hs["z_max"]:.2f} km MSL  '
              f'max={hs["dbz_max"]:.1f} dBZ  '
              f'vox={hs["voxels"]}')

    results.append(scan_result)

print(f'\nDone: {len(results)} scans processed')