import numpy as np
import s3fs
import xarray as xr
from pyproj import Proj, Transformer
from datetime import datetime, timezone

# -- 1. Convert centroids to lat/lon ---------------------------------------
RADAR_LAT =  39.7541
RADAR_LON = -119.4608

# aeqd projection centred on KRGX
aeqd = Proj(proj='aeqd', lat_0=RADAR_LAT, lon_0=RADAR_LON, units='m')
wgs84 = Proj(proj='latlong', datum='WGS84')
to_ll = Transformer.from_proj(aeqd, wgs84, always_xy=True)

# Centroids in km -> metres
plumes = {
    0: dict(x_km=-32.5, y_km=-1.5,  z_km_msl=3.07),
    1: dict(x_km=-24.0, y_km= 4.5,  z_km_msl=2.91),
}
for pid, p in plumes.items():
    lon, lat = to_ll.transform(p['x_km']*1e3, p['y_km']*1e3)
    p['lat'] = lat; p['lon'] = lon
    print(f"Plume {pid}: lat={lat:.4f}  lon={lon:.4f}  "
          f"z={p['z_km_msl']:.2f} km MSL")

# -- 2. Find nearest GOES-18 FDCC file to scan time ------------------------
SCAN_TIME = datetime(2026, 8, 10, 0, 8, 0, tzinfo=timezone.utc)
fs = s3fs.S3FileSystem(anon=True)

# GOES-18 FDCC path: noaa-goes18/ABI-L2-FDCC/YYYY/DOY/HH/
doy = SCAN_TIME.timetuple().tm_yday
hh  = SCAN_TIME.hour
prefix = f"noaa-goes18/ABI-L2-FDCC/{SCAN_TIME.year}/{doy:03d}/{hh:02d}/"
files  = fs.ls(prefix)
print(f"\nGOES FDCC files near {SCAN_TIME.strftime('%H:%MZ')}:")
for f in files:
    print(" ", f.split('/')[-1])

# Pick closest file by start time encoded in filename
# Filename format: OR_ABI-L2-FDCC-M6_G18_s{start}_e{end}_c{created}.nc
def goes_fname_time(path):
    s = path.split('_s')[1][:13]   # YYYYDOYHHMMss
    return datetime.strptime(s[:11], '%Y%j%H%M').replace(tzinfo=timezone.utc)

best = min(files, key=lambda f: abs((goes_fname_time(f)-SCAN_TIME).total_seconds()))
print(f"\nUsing: {best.split('/')[-1]}")

# -- 3. Load and check fire mask at plume locations ------------------------
with fs.open(best) as fobj:
    ds = xr.open_dataset(fobj)
    ds.load()   # force read while the fsspec file handle is still open

# GOES ABI fixed grid -> lat/lon
goes_proj = Proj(proj='geos', lon_0=-137.0, h=35786023.0,
                 sweep='x', ellps='GRS80')
x_rad = ds['x'].values   # radians
y_rad = ds['y'].values
sat_h = 35786023.0
X, Y  = np.meshgrid(x_rad * sat_h, y_rad * sat_h)
to_goes = Transformer.from_proj(wgs84, goes_proj, always_xy=True)

mask  = ds['Mask'].values    # (y, x)
power = ds['Power'].values   # MW
temp  = ds['Temp'].values    # K

print("\nGOES fire mask values: 10-15=fire, 30-35=temporally-filtered fire")
SEARCH_DEG = 0.15   # search radius in degrees (~15 km)

for pid, p in plumes.items():
    # find nearest GOES pixel
    goes_x, goes_y = to_goes.transform(p['lon'], p['lat'])
    dist2 = (X - goes_x)**2 + (Y - goes_y)**2
    iy, ix = np.unravel_index(dist2.argmin(), dist2.shape)

    # check small neighbourhood
    r = 3   # pixel radius
    iy0,iy1 = max(0,iy-r), min(mask.shape[0],iy+r+1)
    ix0,ix1 = max(0,ix-r), min(mask.shape[1],ix+r+1)

    patch_mask  = mask[iy0:iy1, ix0:ix1]
    patch_power = power[iy0:iy1, ix0:ix1]
    patch_temp  = temp[iy0:iy1, ix0:ix1]

    fire_pixels = np.isin(patch_mask, [10,11,12,13,14,15,30,31,32,33,34,35])

    print(f"\nPlume {pid} (lat={p['lat']:.4f}, lon={p['lon']:.4f}):")
    print(f"  Nearest GOES pixel mask value: {patch_mask[r,r]}")
    print(f"  Fire pixels in {2*r+1}x{2*r+1} neighbourhood: "
          f"{fire_pixels.sum()}")
    if fire_pixels.sum() > 0:
        print(f"  Max fire power: {np.nanmax(patch_power[fire_pixels]):.1f} MW")
        print(f"  Max fire temp:  {np.nanmax(patch_temp[fire_pixels]):.1f} K")
    else:
        print("  No fire detection in neighbourhood")

ds.close()
