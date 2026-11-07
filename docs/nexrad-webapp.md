# NEXRAD Webapp and Column-Labeling ML Pipeline

Goal: label pyroconvective "column" structures in 3D NEXRAD reflectivity volumes to train a model, then scale to ~800 fires. Dev case: Camp Fire, KBBX, 2018-11-08 (frames 15:11, 15:21, 15:46 UTC).

## Components
- Radar slicer app: self-contained HTML (Three.js) with 2D slice panels + 3D isosurface panel
  - Discrete threshold level buttons: 20/30/35/40/42/45/47 dBZ (bottom bar)
  - One Three.js mesh per connected component (scipy `label` + skimage `marching_cubes` on isolated component masks)
  - Click a component to label it: Column / Falling column / Not a column; export labels as CSV
  - "AWS S3" side panel: list scans for a time window, select, stream into the app via the server
- `radar_server.py`: Flask backend; fetches from `unidata-nexrad-level2`, grids with pyart Barnes2, returns frame JSON
- `generate_radar_frames.py`: offline frame generator for the app
- `nexrad_grid_loader.py`: S3 fetch (unsigned), gzip-aware, drops sweeps below `MIN_ELEV_DEG=0.6`, Barnes2 gridding, transposes pyart `[z,y,x]` to `[x,y,z]`; returns `(dbz_grid, x_km, y_km, z_km)`
- `smoke_test_real_scan.py`: sanity check, e.g. `python smoke_test_real_scan.py --site KBBX --date 2018-11-08 --hour 15`; prints kept/dropped sweeps and echo centroid to catch axis-transpose bugs
- `threshold_slider_labeler.py`: Dash labeler; `USE_REAL_NEXRAD_DATA`, `REAL_SCAN_SITE`, `REAL_SCAN_KEY` toggles; QGIS-style layer panel (frame > threshold > component eye toggles, parent hides child, state in a `dcc.Store` keyed by path strings); ghost meshes at 10% opacity for other thresholds
- `labels_csv_to_patches.py` (Phase 3, next step): extract fixed-size 3D voxel patches around labeled component centroids; `REAL_SCAN_LOOKUP` maps `(fire_id, timestamp)` to `(site, S3_key)`; synthetic fallback when empty

## Environment
- atar conda env: `/data/Raymond/miniconda3/envs/pyart_312` (Python 3.12, pyart)

## Gotchas
- Marching cubes in JS must use Bourke corner ordering via an explicit `MC_CORNERS` table: `[(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]`; bit-pattern ordering gives a blank 3D panel
- `putImageData` ignores canvas transforms; draw to an offscreen canvas then `drawImage` for zoom/pan
- Barnes2 gridding ~10 s/scan, dominated by S3 download
- Check Barnes2 params (`map_gates_to_grid` vs `map_to_grid`, `h_factor`, `kappa_star`, `gamma`) first if gridded output is all NaN or misplaced
- 16 MB limit applies only to claude.ai-hosted artifacts; the local HTML (~11 MB) is fine
- Be polite to the Unidata S3 bucket (public dataset, avoid hammering it)
- Incremental JS patches have repeatedly blanked the app via runtime errors; check the browser console first
