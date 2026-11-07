# CLAUDE.md: NEXRAD Radar Plume and Boundary Tracking

Context for Claude Code agents working on this project. Owner: Raymond Harding, UNR Wildfire Science and Engineering Lab (advisor: Dr. Neil Lareau).

## Goal
Use NEXRAD Level-II radar (plus VIIRS/MODIS and NIFC IR perimeters) to track wildfire smoke plumes and fire boundaries.

## Sub-projects
- `docs/fireline-pipeline.md`: 2D fireline detection/tracking (repo: github.com/RaymndH/Radar_Firelines)
- `docs/plume-tracking-3d.md`: 3D plume tracking with tobac + custom linker + persistent homology (Camp Fire, KBBX)
- `docs/nexrad-webapp.md`: HTML/Three.js radar slicer + Flask server + Dash labeler for labeling pyroconvective columns (ML training data, target ~800 fires)
- `docs/bug-fire-krgx.md`: Bug Fire 2026 KRGX dual-pol smoke mask, GOES-anchored occupancy tracking, open column-tracking problem

## Code in this repo
- `webapp/`: slicer HTML + sample labels; Python server/loader/labeler scripts are missing, rebuild from `docs/nexrad-webapp.md` (see `webapp/README.md`)
- `bug_fire_krgx/`: Bug Fire scripts, `TASKS.md` (full task log 01 to 25), its own `CLAUDE.md` with gridding params, results JSON/NPZ, final perimeter shapefile; large `grids_*.npz` and frame PNGs stay on the Windows machine at `C:\Users\raymondharding\Research\Runs\krgx_hourly\`

## Next steps
- Webapp pipeline: Phase 3 `labels_csv_to_patches.py` (3D voxel patches around labeled centroids), then scale to ~800 fires
- Bug Fire: Task 26, GOES-seeded marker-controlled watershed for column detection

## Working style
- Debug iteratively from actual error output, not speculative fixes
- Flag deviations from spec, circular logic, and physically implausible results before building on them

## Environment
- Dev: Windows + Anaconda, env `geo_env` (Python 3.11)
- HPC: `atar` (~256 cores, NVIDIA GPUs) and `atar2`; large data on atar, Box for sync, Google Drive for archives
- Key libs: pyart, wradlib, tobac, gudhi, scikit-image, plotly/kaleido, geopandas, GDAL, netCDF4, xarray, boto3, alphashape, shapely 2.x

## Data sources
- NEXRAD: public S3 bucket `unidata-nexrad-level2`; use an unsigned boto3 client (`Config(signature_version=UNSIGNED)`)
- FIRMS/VIIRS: API key read from env var `FIRMS_API_KEY` (never commit keys)
- NIFC IR perimeters for validation

## Conventions
- Projection: radar-centered aeqd; recreate the proj object per radar site (stale proj gives wrong XY when switching fires)
- Gridding: cKDTree-based (~23x faster than `pyart.map.grid_from_radars`)
- Stage 1 output: NetCDF, 30x30 km Cartesian, 100 m resolution, zlib-compressed

## Fires processed
Camp 2018 (KBBX), Bear 2020, Beaver 2014 (KMAX), Lolo Peak 2017 (KMSX), Riverside 2020 (KRTX)

## Known gotchas (read before debugging)
- `.gz` S3 keys: `fetch_bytes()` in `async_preprocess.py` must `gzip.decompress` when key ends in `.gz`, or pyart fails with "unknown compression record"
- `firelines_alg.py` once set `PYTHONBREAKPOINT=0` at import, disabling `breakpoint()`; use `ipdb.set_trace()`
- pdb: list comprehensions fail inside pdb; assign intermediates or use `exec()`
- `alphashape.optimizealpha()` breaks on shapely 2.x (MultiPoint not iterable); use fixed alpha (~50), >=4 point guard, try/except fallback; seed with IR perimeter boundary and keep the largest polygon
- gudhi `min_persistence` does not reliably filter `cofaces_of_persistence_pairs()`; filter births/deaths/persistence manually
- Downsample with `block_reduce(..., np.max)` not striding, to keep peak dBZ
- `marching_cubes` needs each bbox dim > step_size; clamp per region
- Colab: pin `kaleido==0.2.1` before any other install/import of plotly
- Plotly cost: `go.Volume` slow at export, `go.Isosurface` slow interactively, precomputed `Mesh3d` cheapest
