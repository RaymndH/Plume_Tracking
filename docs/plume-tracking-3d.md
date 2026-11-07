# 3D Plume Tracking (tobac)

Case: Camp Fire 2018, KBBX. Notebook: `3d_tobac.ipynb` (Colab, later VS Code + geo_env).
Working dir: `C:\Users\raymondharding\Box\Raymond_General_Materials_2026\Radar_Data_and_Scripts\Runs\2018_Camp_CA\`

## Pipeline
1. `list_scans` / `fetch_radar`: pull scans from S3
2. `grid_scan`: Cartesian gridding with elevation filtering (`min_elev=0.6`)
3. tobac-style 3D feature detection + segmentation
4. `link_irregular_time_predictive`: custom linker, Hungarian assignment (`linear_sum_assignment`) with constant-velocity prediction; replaces `linking_trackpy`, which assumes a fixed dt and breaks on NEXRAD scan gaps / VCP changes
5. Filtering:
   - `filter_features_by_core`: relative dBZ threshold + near-surface core base
   - `is_likely_column`: voxel size, height, persistence, max speed, distance to known fire; separates rooted pyroconvective columns from noise / downwind smoke
6. Topological filtering: `gudhi.CubicalComplex` 0-dim persistent homology on 2D slices (~z=500 m) to pick significant cores; filter pairs manually (see gotchas in CLAUDE.md); downsample with `block_reduce(np.max)`

## Visualization
plotly `go.Volume`, `go.Mesh3d` (from `skimage.measure.marching_cubes`), `go.Surface`, `go.Scatter3d`, `bbox_wireframe_trace`; kaleido for PNG, imageio for MP4
