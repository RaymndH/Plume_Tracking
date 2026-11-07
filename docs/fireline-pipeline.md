# Fireline Pipeline (Radar_Firelines)

Repo: https://github.com/RaymndH/Radar_Firelines

## Architecture
- Stage 1 (`nexrad_preprocessor_optimized.py`, `async_preprocess.py`): NEXRAD Level-II to NetCDF (30x30 km, 100 m, zlib), cKDTree gridding
- Stage 2 (`firelines_alg.py`): fireline tracking with `RegularGridInterpolator`, IoU validation vs VIIRS and NIFC IR perimeters
- `viirs_utils.py`: VIIRS fetch + DBSCAN cumulative perimeters in radar XY (aeqd)
- `fire_pipeline.py`: driver with run types `Default_start`, `VIIRS_start` (VIIRS first detection as start), `IR_start`

## Local paths (Windows)
- Scripts: `C:\Users\raymondharding\Box\Raymond_General_Materials_2026\Radar_Data_and_Scripts\Code\New_Python_Scripts\`
- NetCDF output: `C:\Users\raymondharding\radar_output\`

## Fixes worth remembering
- Beaver Fire (KMAX, 2014-08-02): 125/125 scans failed with "unknown compression record". Objects are outer-gzipped AR2V0006. `stream_scan()` decompressed `.gz` but `fetch_bytes()` did not; add `if key.endswith('.gz'): data = gzip.decompress(data)` to `fetch_bytes()`
- Beaver Fire gif radar perimeter: alphashape fragmented into many polygons; fixed by seeding with IR exterior coords and keeping the max-area polygon
