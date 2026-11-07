import numpy as np
import s3fs, xarray as xr
from pyproj import Proj, Transformer

RADAR_LAT =  39.7541
RADAR_LON = -119.4608

# Same GOES file as Task 07
GOES_FILE = ('noaa-goes18/ABI-L2-FDCC/2026/222/00/'
             'OR_ABI-L2-FDCC-M6_G18_s20262220006181_'
             'e20262220008554_c20262220009040.nc')

fs = s3fs.S3FileSystem(anon=True)
with fs.open(GOES_FILE) as fobj:
    ds = xr.open_dataset(fobj, engine='h5netcdf')
    ds.load()

# GOES ABI fixed grid -> lat/lon
sat_h     = 35786023.0
goes_proj = Proj(proj='geos', lon_0=-137.0, h=sat_h,
                 sweep='x', ellps='GRS80')
wgs84     = Proj(proj='latlong', datum='WGS84')
aeqd      = Proj(proj='aeqd', lat_0=RADAR_LAT, lon_0=RADAR_LON, units='m')

to_ll     = Transformer.from_proj(goes_proj, wgs84, always_xy=True)
to_aeqd   = Transformer.from_proj(wgs84, aeqd, always_xy=True)
aeqd_to_ll = Transformer.from_proj(aeqd, wgs84, always_xy=True)

x_rad = ds['x'].values
y_rad = ds['y'].values
X, Y  = np.meshgrid(x_rad * sat_h, y_rad * sat_h)

mask  = ds['Mask'].values.astype(float)
power = ds['Power'].values.astype(float)
ds.close()

# Fire pixels: mask codes 10-15 (fire) and 30-35 (temporally filtered)
fire_codes = set(range(10,16)) | set(range(30,36))
is_fire = np.isin(mask, list(fire_codes))
print(f"Total GOES fire pixels in scene: {is_fire.sum()}")

# Convert all fire pixel locations to aeqd km
fire_x_goes = X[is_fire]
fire_y_goes = Y[is_fire]
fire_power  = np.where(np.isfinite(power[is_fire]), power[is_fire], 0.0)

fire_lon, fire_lat = to_ll.transform(fire_x_goes, fire_y_goes)
fire_x_km, fire_y_km = [v/1e3 for v in to_aeqd.transform(fire_lon, fire_lat)]

# GOES FDCC covers the whole CONUS scene, so `is_fire` includes every
# wildfire detection currently burning anywhere in the west, not just near
# KRGX. Restrict to the radar analysis domain (krgx_plume.py X_MIN/X_MAX,
# Y_MIN/Y_MAX, in km) before clustering, otherwise the "two sources" split
# below just bisects fires hundreds of km apart.
X_MIN_KM, X_MAX_KM = -50.0, 10.0
Y_MIN_KM, Y_MAX_KM = -10.0, 25.0
in_domain = ((fire_x_km >= X_MIN_KM) & (fire_x_km <= X_MAX_KM) &
             (fire_y_km >= Y_MIN_KM) & (fire_y_km <= Y_MAX_KM))
print(f"Fire pixels within KRGX domain: {in_domain.sum()} / {is_fire.sum()}")

fire_x_km   = fire_x_km[in_domain]
fire_y_km   = fire_y_km[in_domain]
fire_power  = fire_power[in_domain]

# Cluster fire pixels into two groups using the known approximate separation
# (~9 km between the two sources) — simple distance split around midpoint
# Midpoint of the two plume centroids in aeqd km:
mid_x = (-32.5 + -24.0) / 2   # = -28.25 km
mid_y = ( -1.5 +   4.5) / 2   # =   1.50 km

# Split: fire pixels closer to west source vs east source
# Use bearing from midpoint — pixels west of midpoint -> Source A, east -> B
# Simple x-split at mid_x is sufficient given 8.5 km x-separation
src_a = fire_x_km < mid_x
src_b = ~src_a

print(f"\nSource A fire pixels: {src_a.sum()}")
print(f"Source B fire pixels: {src_b.sum()}")

for label, sel, name in [('A', src_a, 'Bug Fire (main)'),
                          ('B', src_b, 'Spot Fire NE')]:
    if sel.sum() == 0:
        print(f"\nSource {label} ({name}): NO FIRE PIXELS — check split")
        continue
    w = fire_power[sel]
    if w.sum() == 0:
        w = np.ones(sel.sum())   # fallback: unweighted if power all zero
    cx_km = np.average(fire_x_km[sel], weights=w)
    cy_km = np.average(fire_y_km[sel], weights=w)
    clon, clat = aeqd_to_ll.transform(cx_km*1e3, cy_km*1e3)  # aeqd km -> lat/lon
    print(f"\nSource {label} - '{name}':")
    print(f"  N fire pixels: {sel.sum()}")
    print(f"  Power-weighted centroid: x={cx_km:.2f} km  y={cy_km:.2f} km")
    print(f"  Lat/lon: {clat:.4f}N  {clon:.4f}W")
    print(f"  Total power: {fire_power[sel].sum():.0f} MW")
    print(f"  Max power:   {fire_power[sel].max():.0f} MW")
