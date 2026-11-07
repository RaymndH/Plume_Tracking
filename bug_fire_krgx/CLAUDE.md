# KRGX Plume Analysis — Project Context

## Environment
- Machine: Windows, conda env `geo_env`
- Working directory: `C:\Users\raymondharding\Research\Runs`
- Python: 3.11, Jupyter (VSCode) + IPython terminal
- Key packages: pyart, boto3, numpy, scipy, sklearn, skimage, plotly, matplotlib, pyproj, s3fs, xarray

## Radar
- Site: KRGX (Reno, NV)
- Lat/Lon: 39.7541, -119.4608
- Altitude: 2530 m MSL (8299 ft)
- Target event: Bug & Stallion Fire, 2026-08-09 ~23:59Z
- Data source: `unidata-nexrad-level2` (primary), `noaa-nexrad-level2` (fallback), anonymous S3
- Cache dir: `./krgx_cache/` — raw Level-II files cached by filename

## Domain (aeqd metres from radar)
- X: -50,000 to -20,000 m (W of radar)
- Y: -10,000 to +25,000 m
- Z: 0 to 8,000 m MSL ceiling; gates filtered to AGL >= 0
- DXY: 250 m, DZ: 250 m
- ELEV_MIN: -0.5°, ELEV_MAX: 19.5°

## Gridding
- Method: `numpy.histogramdd` directly from raw gates — NOT `grid_from_radars`
- Reason: `grid_from_radars` produces a ~5 km spatial offset vs raw gate coordinates
- Gates stored as AGL km internally; converted to MSL km for all plotting
- `z_km_msl = z_km + RADAR_ALT_M / 1e3`
- Plotly meshgrid: `np.meshgrid(x_km, y_km, z_km_msl, indexing='ij')` → shape (nx, ny, nz)
- dbz_grid shape: (nx, ny, nz) matching histogramdd bin order

## Fields
- Reflectivity: primary field, dBZ
- Spectrum width: secondary field (m/s), added to same histogramdd grid
- Both extracted from all sweeps in ELEV_MIN/MAX range

## Key parameters
- DBZ_THRESHOLD: 5–15 dBZ (binary mask for watershed/skeleton)
- GAUSS_SIGMA: 1.0 voxels (applied before binary threshold)
- PEAK_FRAC: 0.4 (peak must be >= 40% of global max to be kept)
- NEIGHBOR_RADIUS_M: 2000 m
- SADDLE_DROP_THRESH: 3.0 dBZ (below this, merge two watershed regions)
- SW_ISO_MIN: 3.0 m/s, SW_ISO_MAX: 12.0 m/s

## Hotspot detection pipeline
1. `histogramdd` → `dbz_grid` (nx, ny, nz)
2. `gaussian_filter` → `dbz_smooth`
3. `dbz_smooth > DBZ_THRESHOLD` → `binary_3d`
4. `dbz_smooth.max(axis=2)` → `col_max_2d` (nx, ny) column max
5. `peak_local_max(col_max_2d, min_distance=5, threshold_abs=8)` → candidate peaks
6. Filter by: `peak_val >= global_max * PEAK_FRAC` AND `n_neighbors >= MIN_NEIGHBOR_VOX`
7. `watershed(-dbz_smooth, markers, mask=binary_3d)` → `labels_ws`
8. Saddle-based merging: if `min(peak-saddle) <= SADDLE_DROP_THRESH` → merge regions
9. Per-region stats: x/y/z centroid, z_max, dbz_max, voxel count

## Skeleton
- `skimage.morphology.skeletonize(binary_3d)` on smoothed grid
- Color by `dbz_smooth[sv_x, sv_y, sv_z]` (intensity at skeleton voxels)
- SW outline: `binary_erosion` on SW binary → boundary scatter in black

## Terrain
- Source: Open-Elevation API (SRTM), cached to `krgx_cache/terrain.npy`
- Shape must match current (ny, nx) — delete cache if DXY changes
- Range in domain: ~1.29–2.68 km MSL (Sierra Nevada)
- Added as LAST plotly trace to avoid corrupting scene bounds
- 2D plots: `ax.contour` overlay at 250m intervals

## GOES
- Satellite: GOES-18 (West)
- Product: `ABI-L2-FDCC` (Fire Detection, CONUS, ~5 min cadence)
- Bucket: `noaa-goes18`, anonymous S3 via s3fs
- Fire mask values: 10–15 (fire), 30–35 (temporally filtered fire)
- Reprojected from GOES geostationary to radar aeqd using pyproj CRS transform
- Key fields: `Mask`, `Power` (MW), `Temp` (K), `Area` (m²)

## NIFC perimeters
- Service: WFIGS `ABI-L2-CMIPF` → actually `WFIGS_Interagency_Perimeters_YearToDate`
- Field names have `poly_` and `attr_` prefixes (schema changed ~2025)
- Fire name field: `attr_IncidentName`
- Bug Fire IRWIN ID: query via `attr_IncidentName LIKE '%BUG%'`
- Only one perimeter in WFIGS (final ~93k acre shape) — no daily IR progression available
- FTP: `ftp.wildfire.gov/protected/incident_specific_data/great_basin/2026_Incidents/2026_Bug/` — IAP PDFs only, no GIS

## Plotly 3D conventions
- All z in km MSL
- Reflectivity isosurface: `colorscale='Blues'`, opacity 0.5
- SW isosurface: `colorscale=[yellow→orange→orangered]`, flatshading=True
- SW outline: `binary_erosion` boundary as black Scatter3d size=2
- Skeleton: `symbol='diamond'`, colored by dBZ intensity, size=5
- Raw gates: `symbol='circle'`, size=2, opacity=0.3, prioritize high-dBZ points
- Terrain: `go.Surface`, colorscale='earth', added LAST
- Vertical exaggeration: z aspectratio = `(Z_range/Y_range) * 3`
- z_plot_min = 1.0 km MSL (shows terrain); z_plot_max = Z_MAX/1e3

## Files
- `krgx_plume.py` — single scan: full pipeline + 3D plotly
- `krgx_hourly.py` — one hour of scans: per-scan 2D PNGs + results list
- `krgx_3d_viz.ipynb` — original 3D viz notebook (timestep slider)
- `krgx_3d_plume.html` — output 3D interactive plot
- `krgx_hourly/` — per-scan PNG output directory

## Known issues / decisions
- `grid_from_radars` discarded: produces ~5 km offset, 18 dBZ value loss vs raw gates
- Terrain cache is resolution-dependent — delete `terrain.npy` when DXY changes
- Plotly Scatter3d does not support per-point opacity — use RGBA color strings instead
- `go.Isosurface` requires fill value (not NaN) for no-data voxels — use `np.nan` (confirmed working)
- `peak_local_max` with large `min_distance` fails when array is thin in z — use 2D column-max projection instead
- Windows Search does not index `Runs` folder by default — use File Explorer search or `dir /s` in terminal
- KRGX at 2530 m MSL — gates below radar site altitude are clutter; filter with `gz_agl >= 0`
## Environment — shell workaround
Claude Code's PowerShell tool runs with `-NoProfile -NonInteractive`, so
`conda activate` never works in tool calls regardless of `conda init` state.
All Claude-Code-run scripts must invoke the interpreter directly:
`AppData\Local\miniforge3\envs\geo_env\python.exe` with
`Library\bin`, `Library\mingw-w64\bin`, and `Scripts` prepended to `$env:PATH`.
This is permanent — do not add "fix conda activate" to any future task.