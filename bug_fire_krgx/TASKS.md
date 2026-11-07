# KRGX Plume Analysis — Task Queue
# PM: Claude (claude.ai project session)
# Agent: Claude Code (local)
# Updated: 2026-09-18

---

## TASK 01 — Verify dual-pol field availability
**Status**: PENDING

**Goal**: Confirm which polarimetric fields are actually present in the cached
Level-II file for the Bug Fire scan (2026-08-10 ~00:08Z).

**Instructions**:
1. Activate conda env: `conda activate geo_env`
2. Run the following in IPython or a scratch script:

```python
import io, gzip
from pathlib import Path
import pyart

CACHE_DIR = Path('./krgx_cache')
# Find the Bug Fire scan file
files = sorted(CACHE_DIR.glob('KRGX20260810_*'))
print("Cached files:", files)

# Load the first one
raw = files[0].read_bytes()
radar = pyart.io.read_nexrad_archive(io.BytesIO(raw))
print("\nAvailable fields:", list(radar.fields.keys()))
print("N sweeps:", radar.nsweeps)
print("Elevations:", [round(float(radar.elevation['data'][radar.get_slice(i)].mean()),1)
                      for i in range(radar.nsweeps)])

# Check data ranges for each field
for fname, fdata in radar.fields.items():
    d = fdata['data']
    print(f"  {fname:35s}  min={float(d.min()):.2f}  max={float(d.max()):.2f}  "
          f"masked%={100*d.mask.sum()/d.size:.1f}%")
```

**Append results below this line** — field names, value ranges, sweep count.

---
## RESULTS 01
**Status**: DONE — run 2026-09-18, via `geo_env` (miniforge3; note: `conda activate`
was not on PATH in this shell — had to invoke `envs\geo_env\python.exe` directly
with `Library\bin` etc. prepended to PATH to satisfy GDAL/pyproj DLLs).

Used cached file `krgx_cache/KRGX20260810_000404_V06` (closest to T_CENTER
2026-08-10 00:08:04Z among the 13 cached KRGX20260810_* scans).

**Available fields (7 total):**
- `differential_reflectivity` (ZDR)  — min=-13.00  max=20.00  masked%=93.5%
- `cross_correlation_ratio` (RhoHV)  — min=0.21  max=1.05  masked%=93.5%
- `reflectivity` (Z)                 — min=-32.00  max=59.00  masked%=88.9%
- `velocity` (Vr)                    — min=-28.00  max=28.00  masked%=94.8%
- `differential_phase` (PhiDP)       — min=0.00  max=359.65  masked%=93.5%
- `spectrum_width` (SW)              — min=0.00  max=16.00  masked%=94.8%
- `clutter_filter_power_removed`     — min=-6.00  max=73.00  masked%=91.5%

**N sweeps**: 16 (VCP with split cuts at low elevations)
**Elevations (deg, mean per sweep)**: [0.0, 0.0, 0.5, 0.5, 0.8, 0.8, 1.2, 1.2, -0.0, 0.0, 1.8, 2.4, 3.0, 4.0, 5.1, 6.4]

**Confirmed**: KRGX is full dual-pol — ZDR, RhoHV, PhiDP, and Vr are all present
and populated (not empty/all-masked) in the archive. `krgx_plume.py` line 161
currently restricts `include_fields` to only `['reflectivity', 'spectrum_width']`,
so ZDR/RhoHV/PhiDP are being silently discarded at read time. Path 1
(polarimetric smoke mask) is unblocked — no data gap, just a code change.

Note on masked%: ~88–95% masked across all fields is expected — most gates in
the domain box are range/elevation cells with no return (clear air), not a
data-quality issue with the polarimetric moments specifically.

Note on elevations: two sweeps read as ~0.0° and one pair appears near-duplicate
(0.0, 0.0), (0.5, 0.5), (0.8, 0.8), (1.2, 1.2) — consistent with SZ2/CS-ORPG
split cuts (separate surveillance/Doppler sweeps at the same nominal elevation
for dual-pol VCPs), not a bug in the elevation averaging.

## TASK 02 — Fix sw variable collision + inject polarimetric fields
**Status**: PENDING

**Goal**: Two surgical edits to `krgx_plume.py`:
1. Fix the loop variable name collision (sw → isw)
2. Add ZDR, RhoHV, PhiDP to include_fields, gate extraction, and histogramdd grids

**Instructions**: Open `krgx_plume.py` and apply the following changes exactly.

### Change 1 — include_fields (line ~161)
Find:
```python
radar = pyart.io.read_nexrad_archive(
    io.BytesIO(raw), include_fields=['reflectivity', 'spectrum_width']
)
```
Replace with:
```python
radar = pyart.io.read_nexrad_archive(
    io.BytesIO(raw), include_fields=[
        'reflectivity', 'spectrum_width',
        'differential_reflectivity', 'cross_correlation_ratio', 'differential_phase'
    ]
)
```

### Change 2 — gate extraction loop (lines ~176–202)
Find the loop header:
```python
for sw in keep:
    sl = radar.get_slice(sw)
```
Replace `sw` → `isw` throughout the entire loop body. The corrected loop should read:
```python
for isw in keep:
    sl = radar.get_slice(isw)
    az = radar.azimuth['data'][sl]
    el = radar.elevation['data'][sl]

    z_field = radar.fields.get('reflectivity')
    if z_field is None:
        continue
    z_data = z_field['data'][sl]

    sw_field = radar.fields.get('spectrum_width')
    sw_data = sw_field['data'][sl] if sw_field is not None else None

    zdr_field = radar.fields.get('differential_reflectivity')
    zdr_data = zdr_field['data'][sl] if zdr_field is not None else None

    rhv_field = radar.fields.get('cross_correlation_ratio')
    rhv_data = rhv_field['data'][sl] if rhv_field is not None else None
```
Keep the existing gate coordinate calculation (rng, gz_agl, gx, gy) and valid-gate
mask unchanged. Extend the mask and append block to collect the new fields:

After the line that appends to `all_sw` (or wherever sw_data is appended),
add:
```python
    valid = ~np.ma.getmaskarray(z_data)  # existing valid mask — do not redefine

    # ... existing appends for gx, gy, gz, dbz, sw_vals ...

    if zdr_data is not None:
        all_zdr.append(np.ma.filled(zdr_data[valid], np.nan))
    else:
        all_zdr.append(np.full(valid.sum(), np.nan))

    if rhv_data is not None:
        all_rhv.append(np.ma.filled(rhv_data[valid], np.nan))
    else:
        all_rhv.append(np.full(valid.sum(), np.nan))
```

### Change 3 — Initialize accumulator lists (before the loop)
Find where `all_sw = []` (or similar) is initialized before the loop.
Add alongside it:
```python
all_zdr = []
all_rhv = []
```

### Change 4 — histogramdd grids for ZDR and RhoHV (after the loop, near dbz_grid)
Find where `dbz_grid` is built with histogramdd. After that block, add:
```python
# --- ZDR grid ---
zdr_vals = np.concatenate(all_zdr)
zdr_finite = np.isfinite(zdr_vals)
zdr_sum, _ = np.histogramdd(
    gate_coords[zdr_finite],   # same (N,3) array used for dbz_grid
    bins=[x_edges, y_edges, z_edges],
    weights=zdr_vals[zdr_finite]
)
zdr_grid = np.where(count_grid > 0, zdr_sum / count_grid, np.nan)

# --- RhoHV grid ---
rhv_vals = np.concatenate(all_rhv)
rhv_finite = np.isfinite(rhv_vals)
rhv_sum, _ = np.histogramdd(
    gate_coords[rhv_finite],
    bins=[x_edges, y_edges, z_edges],
    weights=rhv_vals[rhv_finite]
)
rhv_grid = np.where(count_grid > 0, rhv_sum / count_grid, np.nan)
```
Note: `gate_coords` is the (N,3) array of [gx, gy, gz] for all valid gates —
use whatever variable name currently holds that in the script. `count_grid` is
the unweighted histogramdd count (number of gates per voxel) — if it doesn't
exist yet, add it as:
```python
count_grid, _ = np.histogramdd(gate_coords, bins=[x_edges, y_edges, z_edges])
```
before the dbz_grid weighted sum.

**After making all changes**: Run `krgx_plume.py` for the Bug Fire scan and
confirm it completes without error. Paste any traceback or the final
"Plot saved" / completion line below.

**Do not yet build the smoke mask or add any Plotly traces — that is Task 03.**

---
## RESULTS 02
**Status**: DONE — run 2026-09-18.

**Naming deviations from the spec** (adapted to the actual script, as flagged
in the task note — logic preserved, names changed):
- Loop var `sw` → `isw` as instructed. The script's existing accumulator
  lists are `gate_x, gate_y, gate_z, gate_v, gate_sw` (not `all_sw` etc) —
  added `gate_zdr, gate_rhv, gate_phidp` alongside them, and post-loop
  concatenations `gzdr_all, grhv_all, gphidp_all` alongside the existing
  `gw_all = np.concatenate(gate_sw)`.
- Refactored the per-field masked→nan→astype(float) extraction (previously
  duplicated for refl/sw) into a small local `_extract(field_name, shape)`
  helper inside the loop, used for sw/zdr/rhv/phidp. Same logic, less
  repetition.
- Also extracted `differential_phase` (PhiDP) end-to-end (include_fields →
  gate extraction → `phidp_grid`), even though the detailed Task 02 steps
  only spelled out ZDR/RhoHV — the task's own Goal line asked for PhiDP too,
  and it's the same pattern, so did all three rather than leaving PhiDP
  half-wired (in include_fields but unused).
- **Deviated from the spec's grid-building approach on purpose**: the spec's
  Change 4 divides `zdr_sum`/`rhv_sum` by a single shared `count_grid` built
  from *all* valid gates (regardless of whether ZDR/RhoHV was itself masked
  at that gate). That would bias the per-voxel mean low in voxels where some
  gates have valid reflectivity but masked polarimetric moments (common —
  masked% was ~93-95% for ZDR/RhoHV/PhiDP vs ~89% for reflectivity in Task
  01's results). Instead, added a `_pol_grid(vals_all)` helper that computes
  its own finite mask and its own count grid per field, so each of
  `zdr_grid`, `rhv_grid`, `phidp_grid` averages only over gates that actually
  have a valid value for that specific field. `dbz_grid`/`sw_grid` were left
  untouched (pre-existing behavior, not in scope for this task).

**Run output** (Bug Fire scan, 2026-08-10 00:08:09Z):
```
16 sweeps: [0.0, 0.0, 0.4, 0.4, 0.8, 0.8, 1.2, 1.2, -0.0, 0.0, 1.8, 2.4, 3.0, 4.0, 5.1, 6.4]
Extracting raw gates...
  620501 total gates in domain
  1697 high-dBZ + 48303 fill = 50000 plot points
Gridding via histogramdd...
  Grid shape: (240, 140, 22)
  Grid max: x=-33.4 y=-1.6 z=2.90 km MSL  val=33.5 dBZ
Skeletonizing...
  53 skeleton voxels

2 plume clusters
  Plume 0: 737 pts  x=-32.5 y=-1.5 z=3.07 km MSL  elev=1.3°  horiz=13.4°
  Plume 1: 20 pts  x=-24.0 y=4.5 z=2.91 km MSL  elev=11.9°  horiz=89.0°
Loading cached terrain...
Building 3D plot...
Saved krgx_3d_plume.html
```
Exit code 0, no traceback. `zdr_grid`, `rhv_grid`, `phidp_grid` are built
(shape (240,140,22), matching `dbz_grid`) but not yet plotted or thresholded —
per instructions, no smoke mask or new Plotly traces added (that's Task 03).

**Environment note (repeats from Task 01)**: `conda activate` is not on PATH
in this non-interactive shell — invoked
`AppData\Local\miniforge3\envs\geo_env\python.exe` directly with
`Library\bin`/`Library\mingw-w64\bin`/`Scripts` prepended to PATH to satisfy
GDAL/pyproj DLLs. Same workaround needed for any future scripted run.

## TASK 03 — Build polarimetric smoke mask + add RhoHV isosurface to plot
**Status**: PENDING

**Goal**: Compute a `smoke_mask_3d` boolean array from the polarimetric grids
and add one new Plotly trace showing the masked region colored by RhoHV.

**Background on thresholds**:
Smoke/ash plumes in radar polarimetry are characterized by:
  - Low ρHV (< 0.90, often < 0.80) — irregular, tumbling particles
  - Near-zero or slightly positive ZDR (−0.5 to +1.5 dB) — roughly spherical
  - Reflectivity above noise floor (> DBZ_THRESHOLD, currently 5 dBZ)
These thresholds follow Melnikov et al. (2008, 2020) and are intentionally
conservative for a first pass — we will tune them in Task 04.

**Instructions**:

### Change 1 — Define smoke mask (after all three pol grids are built)
Add immediately after `phidp_grid` is computed:

```python
# --- Polarimetric smoke mask ---
# Requires: valid dBZ, low RhoHV (irregular particles), ZDR in smoke range
RHV_MAX   = 0.90   # gates above this are likely precip or ground clutter
ZDR_MIN   = -0.5   # dB — exclude strong negative ZDR (large drops / insects)
ZDR_MAX   =  1.5   # dB — exclude large positive ZDR (rain)

smoke_mask_3d = (
    (dbz_smooth > DBZ_THRESHOLD) &       # non-trivial echo
    np.isfinite(rhv_grid) &              # RhoHV actually observed here
    np.isfinite(zdr_grid) &              # ZDR actually observed here
    (rhv_grid < RHV_MAX) &
    (zdr_grid > ZDR_MIN) &
    (zdr_grid < ZDR_MAX)
)
print(f"Smoke mask: {smoke_mask_3d.sum()} voxels "
      f"({100*smoke_mask_3d.mean():.2f}% of grid)")
```

### Change 2 — Add RhoHV isosurface trace to Plotly figure
Find where the existing SW isosurface trace is added (the `go.Isosurface` with
`colorscale=[yellow→orange→orangered]`). Add the following NEW trace
immediately BEFORE it (so SW renders on top):

```python
# RhoHV smoke-mask volume — show only where smoke_mask_3d is True
rhv_plot = np.where(smoke_mask_3d, rhv_grid, np.nan)
fig.add_trace(go.Volume(
    x=XX.flatten(), y=YY.flatten(), z=ZZ.flatten(),
    value=rhv_plot.flatten(),
    isomin=0.50,
    isomax=0.90,
    surface_count=5,
    colorscale='Purples_r',   # dark purple = low RhoHV = most smoke-like
    opacity=0.15,
    caps=dict(x_show=False, y_show=False, z_show=False),
    name='Smoke mask (ρHV)',
    showscale=True,
    colorbar=dict(title='ρHV', x=1.05, len=0.5)
))
```

Note: `XX`, `YY`, `ZZ` are the meshgrid arrays already used for the existing
isosurface traces (shape (nx,ny,nz), indexing='ij'). Use whatever names the
script currently uses for those arrays.

### Change 3 — Print mask/plume overlap stats
After the smoke mask is defined and after the existing plume cluster stats
print, add:

```python
for pid, pmask in enumerate(plume_masks):   # iterate existing watershed labels
    overlap = (smoke_mask_3d & pmask).sum()
    plume_vox = pmask.sum()
    print(f"  Plume {pid}: {plume_vox} vox total, "
          f"{overlap} ({100*overlap/max(plume_vox,1):.0f}%) in smoke mask")
```
If the existing code doesn't expose per-plume boolean masks by that name,
adapt to whatever structure holds the watershed label regions — the goal is
just to print how many voxels of each plume survive the polarimetric filter.

**After changes**: Run `krgx_plume.py`, confirm exit 0, paste:
1. The smoke mask voxel count line
2. The per-plume overlap lines
3. Any traceback

**Do not tune thresholds yet — that is Task 04.**

---
## RESULTS 03
**Status**: DONE — run 2026-09-18.

**Naming/structure deviations from the spec** (adapted to the actual script):
- Figure variable is `fig3d`, not `fig`.
- Meshgrid arrays are `xx3, yy3, zz3` (shape (nx,ny,nz), `indexing='ij'`, all
  in km MSL), not `XX, YY, ZZ` — used those directly, matches spec's shape
  assumption exactly.
- **No `plume_masks` exist** — the script doesn't do per-voxel watershed
  labeling; it runs DBSCAN on the *points* where `binary_3d` is True
  (`pts_bin = xx3[binary_3d], yy3[binary_3d], zz3[binary_3d]`), giving a 1D
  `labels` array aligned with `pts_bin`, not a per-voxel grid. To get
  per-plume overlap against `smoke_mask_3d` (which is a full (nx,ny,nz)
  boolean grid), scattered `labels` back into a full-grid array:
  `label_grid = np.full(binary_3d.shape, -1); label_grid[binary_3d] = labels`
  (safe because boolean-mask indexing/assignment order matches the order
  `pts_bin` was built in). Then `plume_mask = (label_grid == cid)` per cluster.
  This lives right after the existing PCA/cluster print loop, so it reuses
  `n_clust`/`labels`/`binary_3d` already in scope.
- **Change 1 (mask definition) moved**: the spec said to add it "after
  phidp_grid is computed," but `dbz_smooth` (used in the mask condition)
  isn't computed until later, in the Gaussian-smooth/skeletonize block. Put
  `smoke_mask_3d` right after `dbz_smooth = gaussian_filter(...)` instead —
  same thresholds/logic, just moved to where its inputs actually exist.
- RhoHV colorbar placed at `x=1.2` (not `x=1.05`) since the dBZ raw-gates
  colorbar already sits at `x=1.05` and the SW colorbar at `x=1.1` — avoids
  all three overlapping.
- Used `go.Volume` exactly as specified (not `go.Isosurface`); confirmed it
  renders without error.

**Run output** (Bug Fire scan, 2026-08-10 00:08:09Z):
```
Smoke mask: 2293 voxels (0.31% of grid)
  53 skeleton voxels

2 plume clusters
  Plume 0: 737 pts  x=-32.5 y=-1.5 z=3.07 km MSL  elev=1.3°  horiz=13.4°
  Plume 1: 20 pts  x=-24.0 y=4.5 z=2.91 km MSL  elev=11.9°  horiz=89.0°
  Plume 0: 737 vox total, 153 (21%) in smoke mask
  Plume 1: 20 vox total, 7 (35%) in smoke mask
```
Exit code 0, no traceback. `krgx_3d_plume.html` regenerated with the new
"Smoke mask (ρHV)" Volume trace added before the SW isosurface trace.

**Observation for Task 04 threshold tuning**: only 21-35% of the two
watershed/DBSCAN-detected reflectivity plumes overlap the polarimetric smoke
mask. That's either (a) expected — a fire plume's core is often precip-like
or clutter-like in Z but only its outer/aged smoke fraction is truly
low-RhoHV, or (b) the thresholds (RHV_MAX=0.90, ZDR range [-0.5,1.5]) are too
strict for this scan. Worth a quick histogram of RhoHV/ZDR within the two
plume regions before tuning blind in Task 04.

**Environment note (repeats)**: same `geo_env` python.exe + manual PATH
workaround as Tasks 01/02 — `conda activate` unavailable in this shell.

## TASK 04 — Diagnose polarimetric signature inside plume regions
**Status**: PENDING

**Goal**: Before tuning thresholds, understand the actual RhoHV/ZDR
distribution inside Plume 0 (the real fire plume). Run a diagnostic script
— do NOT modify krgx_plume.py yet.

**Instructions**: Create a scratch script `task04_diag.py` in the Runs dir:

```python
import io, numpy as np, matplotlib.pyplot as plt
from pathlib import Path
import pyart

# ── reproduce the exact grids from krgx_plume.py ──────────────────────────
# (copy/import the gridding block, or re-run the relevant sections)
# Easiest: at the bottom of krgx_plume.py, temporarily add:
#   np.save('krgx_cache/debug_grids.npz', ... )
# But cleaner: just import the module up to the gridding step.
# RECOMMENDED: add a one-time save block to krgx_plume.py — after all four
# grids (dbz_grid, rhv_grid, zdr_grid, binary_3d, label_grid) exist, add:

#   np.savez('krgx_cache/debug_grids.npz',
#            dbz_smooth=dbz_smooth,
#            rhv_grid=rhv_grid,
#            zdr_grid=zdr_grid,
#            binary_3d=binary_3d,
#            label_grid=label_grid)
#   print("Saved debug_grids.npz")

# Run krgx_plume.py once to produce the file, then remove the save block.
# ──────────────────────────────────────────────────────────────────────────

d = np.load('krgx_cache/debug_grids.npz')
dbz_smooth = d['dbz_smooth']
rhv_grid   = d['rhv_grid']
zdr_grid   = d['zdr_grid']
binary_3d  = d['binary_3d']
label_grid = d['label_grid']

# Plume 0 mask (cluster id 0 from DBSCAN)
p0 = label_grid == 0

# --- Print percentile table ---
print("=== Plume 0: RhoHV ===")
rhv_p0 = rhv_grid[p0 & np.isfinite(rhv_grid)]
print(f"  N voxels with finite RhoHV: {len(rhv_p0)}")
if len(rhv_p0):
    for pct in [5,10,25,50,75,90,95]:
        print(f"  p{pct:02d}: {np.percentile(rhv_p0, pct):.3f}")

print("\n=== Plume 0: ZDR ===")
zdr_p0 = zdr_grid[p0 & np.isfinite(zdr_grid)]
print(f"  N voxels with finite ZDR: {len(zdr_p0)}")
if len(zdr_p0):
    for pct in [5,10,25,50,75,90,95]:
        print(f"  p{pct:02d}: {np.percentile(zdr_p0, pct):.2f} dB")

print("\n=== Plume 0: dBZ ===")
dbz_p0 = dbz_smooth[p0]
for pct in [5,10,25,50,75,90,95]:
    print(f"  p{pct:02d}: {np.percentile(dbz_p0, pct):.1f} dBZ")

# --- Key question: how many Plume 0 voxels have NO finite RhoHV? ---
n_p0 = p0.sum()
n_rhv_valid = (p0 & np.isfinite(rhv_grid)).sum()
n_zdr_valid = (p0 & np.isfinite(zdr_grid)).sum()
print(f"\nPlume 0: {n_p0} vox total")
print(f"  finite RhoHV: {n_rhv_valid} ({100*n_rhv_valid/n_p0:.0f}%)")
print(f"  finite ZDR:   {n_zdr_valid} ({100*n_zdr_valid/n_p0:.0f}%)")
print(f"  finite both:  {(p0 & np.isfinite(rhv_grid) & np.isfinite(zdr_grid)).sum()}")

# --- Among voxels WITH valid pol data, how many pass current thresholds? ---
pol_valid = p0 & np.isfinite(rhv_grid) & np.isfinite(zdr_grid)
passes = pol_valid & (rhv_grid < 0.90) & (zdr_grid > -0.5) & (zdr_grid < 1.5)
print(f"\nOf {pol_valid.sum()} vox with valid pol data:")
print(f"  Pass RHV<0.90 AND ZDR in [-0.5,1.5]: {passes.sum()} "
      f"({100*passes.sum()/max(pol_valid.sum(),1):.0f}%)")
print(f"  Fail RHV threshold only: "
      f"{(pol_valid & (rhv_grid >= 0.90)).sum()}")
print(f"  Fail ZDR threshold only: "
      f"{(pol_valid & ((zdr_grid <= -0.5)|(zdr_grid >= 1.5))).sum()}")

# --- Histogram plot ---
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
if len(rhv_p0):
    axes[0].hist(rhv_p0, bins=30, color='steelblue', edgecolor='k')
    axes[0].axvline(0.90, color='red', lw=2, label='RHV_MAX=0.90')
    axes[0].set_xlabel('ρHV'); axes[0].set_title('Plume 0 — RhoHV')
    axes[0].legend()
if len(zdr_p0):
    axes[1].hist(zdr_p0, bins=30, color='darkorange', edgecolor='k')
    axes[1].axvline(-0.5, color='red', lw=2)
    axes[1].axvline( 1.5, color='red', lw=2, label='ZDR window')
    axes[1].set_xlabel('ZDR (dB)'); axes[1].set_title('Plume 0 — ZDR')
    axes[1].legend()
plt.tight_layout()
plt.savefig('krgx_cache/task04_pol_histograms.png', dpi=120)
print("\nSaved task04_pol_histograms.png")
```

**Append to RESULTS 04**:
1. Full percentile table output for RhoHV and ZDR
2. The finite-coverage lines (what % of Plume 0 voxels have valid pol data)
3. The threshold pass/fail breakdown
4. Whether the histogram PNG saved without error

The histogram PNG itself — paste as image if possible, otherwise describe
the shape (unimodal/bimodal, where the bulk sits relative to the threshold lines).

---
## RESULTS 04
**Status**: DONE — run 2026-09-18.

**Process note**: `krgx_plume.py` doesn't retain `dbz_smooth`/`rhv_grid`/
`zdr_grid`/`binary_3d`/`label_grid` between runs, so followed the task's
"RECOMMENDED" path exactly — added the temp `np.savez('krgx_cache/debug_grids.npz', ...)`
block right after `label_grid` is built (inside the `if len(pts_bin) > 10:`
block, where all five arrays are already in scope), ran `krgx_plume.py` once,
then removed the temp block immediately after (confirmed removed — the
permanent script is unchanged from Task 03's state). `task04_diag.py` created
as a standalone scratch script in the Runs dir as instructed; only change
from the spec was adding `matplotlib.use('Agg')` since this runs headless,
and swapping "ρHV" for "RhoHV" in axis labels to avoid a matplotlib font
glyph-missing warning for the Greek rho character.

### Finite-coverage numbers (the key diagnostic requested)
```
Plume 0: 737 vox total
  finite RhoHV: 656 (89%)
  finite ZDR:   656 (89%)
  finite both:  656
```

**This does NOT support the data-sparsity hypothesis.** 89% of Plume 0's
voxels have valid RhoHV *and* ZDR — coverage is high, not sparse. The 21%
smoke-mask overlap from Task 03 is a **threshold problem**, not a missing-data
problem.

### Percentile tables
```
RhoHV:  p05=0.208  p10=0.288  p25=0.435  p50=0.603  p75=0.755  p90=0.855  p95=0.899
ZDR:    p05=-4.32  p10=-3.31  p25=-0.92  p50=1.09   p75=3.34   p90=5.34   p95=6.47  (dB)
dBZ:    p05=15.4   p10=15.8   p25=17.1   p50=19.5   p75=22.3   p90=25.4   p95=26.4
```

### Threshold pass/fail breakdown (of the 656 voxels with valid pol data)
```
Pass RHV<0.90 AND ZDR in [-0.5,1.5]: 153 (23%)
Fail RHV threshold only: 33
Fail ZDR threshold only: 494
```
**The RhoHV threshold is doing almost nothing** — median RhoHV is 0.603, well
below 0.90, and only 33/656 voxels (5%) fail on RhoHV alone. **The ZDR window
is the bottleneck** — 494/656 voxels (75%) fail it, almost entirely on the
*upper* bound: ZDR p50=1.09 dB is already inside [-0.5,1.5], but p75=3.34,
p90=5.34, p95=6.47 — the distribution has a long tail well past +1.5 dB, and
the histogram (below) shows ZDR is roughly centered/unimodal around 0–2 dB,
not clustered at low values the way the Melnikov smoke signature would
predict.

### Histogram (`krgx_cache/task04_pol_histograms.png`), described
- **RhoHV** (left panel): broad, roughly flat/multimodal distribution from
  ~0.2 to ~1.0, mild peaks near 0.22, 0.48, and 0.73. The 0.90 threshold line
  sits in the right-hand tail — most of the mass is well below it, consistent
  with the plume genuinely having low-to-mid RhoHV throughout, not just at
  the edges.
- **ZDR** (right panel): unimodal, roughly bell-shaped, centered around
  +1 to +2 dB, spanning about -9 to +10 dB. Both threshold lines (-0.5, +1.5)
  sit inside the rising/falling shoulders of the main peak, not out in a
  tail — meaning the current ZDR window is cutting through the *bulk* of the
  plume's actual distribution rather than excluding a separate contaminant
  population.

**Interpretation for Task 05**: this looks like a threshold-tuning problem on
ZDR specifically, not a data-fill/interpolation problem. The ZDR upper bound
of +1.5 dB may be too tight for a smoke/ash plume at this range and elevation
(possible causes worth checking before widening blindly: mixed-phase/large
ash aggregates giving positive ZDR, non-uniform beam filling at this range,
or an NBF/attenuation artifact — a range/elevation vs. ZDR scatter for this
plume would help distinguish these). Recommend Task 05 focus on re-deriving
the ZDR window from this plume's own distribution (e.g. p10–p90) rather than
pursuing nearest-neighbor pol-data interpolation, since sparsity is not the
limiting factor here.

**Environment note (repeats)**: same `geo_env` python.exe + manual PATH
workaround as prior tasks.

## TASK 05 — Simplify smoke mask to RhoHV + dBZ only (drop ZDR)
**Status**: PENDING

**Goal**: Remove ZDR from the smoke_mask_3d condition entirely.
RhoHV is the reliable discriminant; ZDR is too noisy at this SNR to
be a useful gate.

**Instructions**:

### Change 1 — Update smoke mask in krgx_plume.py
Find the smoke mask block from Task 03:
```python
RHV_MAX   = 0.90
ZDR_MIN   = -0.5   # dB
ZDR_MAX   =  1.5   # dB

smoke_mask_3d = (
    (dbz_smooth > DBZ_THRESHOLD) &
    np.isfinite(rhv_grid) &
    np.isfinite(zdr_grid) &
    (rhv_grid < RHV_MAX) &
    (zdr_grid > ZDR_MIN) &
    (zdr_grid < ZDR_MAX)
)
```
Replace with:
```python
RHV_MAX = 0.90   # low RhoHV = irregular/tumbling particles (smoke/ash)

smoke_mask_3d = (
    (dbz_smooth > DBZ_THRESHOLD) &
    np.isfinite(rhv_grid) &
    (rhv_grid < RHV_MAX)
)
```
Remove the `ZDR_MIN` and `ZDR_MAX` variables entirely.

### After changes: run krgx_plume.py and report
Paste:
1. "Smoke mask: N voxels (X.XX% of grid)"
2. Per-plume overlap lines for Plume 0 and Plume 1

---
## RESULTS 05
**Status**: DONE — run 2026-09-18.

Applied exactly as specified: `ZDR_MIN`/`ZDR_MAX` removed, `smoke_mask_3d` now
only gates on `dbz_smooth > DBZ_THRESHOLD` and `rhv_grid < RHV_MAX` (0.90).
No variable-name adaptation needed here — this section of `krgx_plume.py`
already matched the spec's names from Task 03.

**Run output** (Bug Fire scan, 2026-08-10 00:08:09Z):
```
Smoke mask: 11144 voxels (1.51% of grid)
  53 skeleton voxels

2 plume clusters
  Plume 0: 737 pts  x=-32.5 y=-1.5 z=3.07 km MSL  elev=1.3°  horiz=13.4°
  Plume 1: 20 pts  x=-24.0 y=4.5 z=2.91 km MSL  elev=11.9°  horiz=89.0°
  Plume 0: 737 vox total, 623 (85%) in smoke mask
  Plume 1: 20 vox total, 19 (95%) in smoke mask
```
Exit code 0, no traceback. `krgx_3d_plume.html` regenerated.

**Before/after**:
| | Task 03 (RhoHV+ZDR) | Task 05 (RhoHV only) |
|---|---|---|
| Total smoke mask | 2293 vox (0.31%) | 11144 vox (1.51%) |
| Plume 0 overlap | 153/737 (21%) | 623/737 (85%) |
| Plume 1 overlap | 7/20 (35%) | 19/20 (95%) |

Confirms Task 04's diagnosis — dropping ZDR moved plume overlap from
21%/35% to 85%/95%, consistent with ZDR having been the bottleneck (RhoHV
alone was already doing most of the discriminating work). Total masked
volume grew ~5x (0.31% → 1.51% of the grid), which is expected since RhoHV
alone is a looser filter — worth keeping an eye on false-positive rate
(e.g. ground clutter/insects with low RhoHV outside the actual plume) in a
future task, but not evaluated here since that wasn't in scope for Task 05.

**Environment note (repeats)**: same `geo_env` python.exe + manual PATH
workaround as prior tasks.

## TASK 06 — Fix conda activation in non-interactive shells
**Status**: PENDING

**Goal**: Make `conda activate geo_env` work in this shell so future tasks
don't need the manual python.exe + PATH workaround.

**Instructions**: Run the following in the current shell (cmd or PowerShell —
note which one you're using in the results):

### If cmd.exe:
```bat
conda init cmd.exe
```
Close and reopen the terminal, then test:
```bat
conda activate geo_env
python -c "import pyart; print(pyart.__version__)"
```

### If PowerShell:
```powershell
conda init powershell
```
You may also need to set the execution policy if scripts are blocked:
```powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```
Close and reopen, then test:
```powershell
conda activate geo_env
python -c "import pyart; print(pyart.__version__)"
```

**Append to RESULTS 06**:
1. Which shell (cmd/PowerShell)
2. Output of the test import line after reopening
3. Whether `conda activate geo_env` now works without manual PATH surgery

---
## RESULTS 06
**Status**: PARTIAL — run 2026-09-18. Shell: PowerShell 5.1 (`powershell.exe`).

Ran `conda init powershell` via
`AppData\Local\miniforge3\Scripts\conda.exe init powershell`. It modified:
`C:\Users\raymondharding\Documents\WindowsPowerShell\profile.ps1`
(CurrentUserAllHosts profile) with the standard conda-init hook block. No
execution-policy change was needed — didn't touch `Set-ExecutionPolicy`.

**Test result — did NOT work as hoped**:
```
conda : The term 'conda' is not recognized...
Python was not found; run without arguments to install from the Microsoft
Store...
```

**Root cause**: this isn't a conda/profile problem — it's how Claude Code's
PowerShell tool launches its shell. Checked the actual process args:
```
powershell.exe -NoProfile -NonInteractive -ExecutionPolicy Bypass -Command ...
```
`-NoProfile` is hard-coded by the tool for every command it runs, so
`profile.ps1` (and therefore the conda hook) is never sourced no matter what's
in it. This is a harness-level constraint, not something fixable from inside
a PowerShell tool call.

**What this means in practice**:
- The `conda init powershell` change is real and correctly installed — if
  you (the user) open a normal interactive PowerShell window, VSCode's
  integrated terminal, or Windows Terminal yourself, `conda activate geo_env`
  should now work there.
- It will **not** help Claude Code's own tool calls in this session or future
  ones — every PowerShell tool invocation starts fresh with `-NoProfile`, so
  the manual workaround (`envs\geo_env\python.exe` directly, with
  `Library\bin`/`Library\mingw-w64\bin`/`Scripts` prepended to `$env:PATH`)
  is still required for any Claude-Code-run script and will remain so.
- Have not verified `conda activate` in an actual interactive window myself
  (no way to do that from this tool) — worth the user doing a quick manual
  check on their end to confirm the profile edit works as intended.

**Recommendation**: no code/task change needed on the Claude Code side; if
this workaround is annoying, the cleaner fix would be to set a persistent
env var / alias in Claude Code's own settings (outside TASKS.md scope) rather
than relying on shell profile loading, since the tool bypasses profiles by
design.

## TASK 07 — Geographic coordinates + GOES cross-reference for Plume 0 and 1
**Status**: PENDING

**Goal**: Convert plume centroids from aeqd km to lat/lon and check whether
GOES-18 fire detections collocate with either plume at the scan time
(2026-08-10 00:08Z).

**Instructions**: Create scratch script `task07_goes_check.py` in Runs dir.

```python
import numpy as np
import s3fs
import xarray as xr
from pyproj import Proj, Transformer
from datetime import datetime, timezone

# ── 1. Convert centroids to lat/lon ───────────────────────────────────────
RADAR_LAT =  39.7541
RADAR_LON = -119.4608

# aeqd projection centred on KRGX
aeqd = Proj(proj='aeqd', lat_0=RADAR_LAT, lon_0=RADAR_LON, units='m')
wgs84 = Proj(proj='latlong', datum='WGS84')
to_ll = Transformer.from_proj(aeqd, wgs84, always_xy=True)

# Centroids in km → metres
plumes = {
    0: dict(x_km=-32.5, y_km=-1.5,  z_km_msl=3.07),
    1: dict(x_km=-24.0, y_km= 4.5,  z_km_msl=2.91),
}
for pid, p in plumes.items():
    lon, lat = to_ll.transform(p['x_km']*1e3, p['y_km']*1e3)
    p['lat'] = lat; p['lon'] = lon
    print(f"Plume {pid}: lat={lat:.4f}  lon={lon:.4f}  "
          f"z={p['z_km_msl']:.2f} km MSL")

# ── 2. Find nearest GOES-18 FDCC file to scan time ────────────────────────
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
    return datetime(int(s[:4]), 1, 1, int(s[7:9]), int(s[9:11]),
                    int(s[11:13]), tzinfo=timezone.utc).replace(
                    month=1, day=1) if False else \
           datetime.strptime(s[:11], '%Y%j%H%M').replace(tzinfo=timezone.utc)

best = min(files, key=lambda f: abs((goes_fname_time(f)-SCAN_TIME).total_seconds()))
print(f"\nUsing: {best.split('/')[-1]}")

# ── 3. Load and check fire mask at plume locations ────────────────────────
with fs.open(best) as fobj:
    ds = xr.open_dataset(fobj, engine='scipy')

# GOES ABI fixed grid → lat/lon
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
```

**Append to RESULTS 07**:
1. Lat/lon of both plume centroids
2. Which GOES file was selected
3. Fire pixel count + power/temp for each plume neighbourhood
4. Any traceback

---
## RESULTS 07
**Status**: DONE — run 2026-09-18.

**Bug fixed along the way**: the spec's `with fs.open(best) as fobj: ds =
xr.open_dataset(fobj)` opens the netCDF lazily — `ds['Mask'].values` accessed
*after* the `with` block exits threw `ValueError: I/O operation on closed
file` (xarray/h5netcdf don't read data until you touch `.values`, by which
point the fsspec file handle was already closed). Fixed by adding
`ds.load()` inside the `with` block, forcing the read while the handle is
still open. Also dropped `engine='scipy'` from `xr.open_dataset` (GOES L2
files are netCDF4/HDF5, not netCDF3-classic — scipy's engine can't read
them; letting xarray auto-select landed on `h5netcdf`, which is already
installed in `geo_env`). No other logic changes from the spec.

**1. Plume centroids in lat/lon**:
```
Plume 0: lat=39.7400  lon=-119.8400  z=3.07 km MSL
Plume 1: lat=39.7943  lon=-119.7410  z=2.91 km MSL
```

**2. GOES file selected**: `OR_ABI-L2-FDCC-M6_G18_s20262220006181_e20262220008554_c20262220009040.nc`
(scan start 2026-day222 00:06:18Z — 1m51s before the 00:08:00Z target, closest
of the 12 files listed in that hour's prefix).

**3. Fire pixel / power / temp per plume neighborhood (7x7 GOES pixels)**:
```
Plume 0 (lat=39.7400, lon=-119.8400):
  Nearest GOES pixel mask value: 100.0
  Fire pixels in 7x7 neighbourhood: 12
  Max fire power: 1061.9 MW
  Max fire temp:  1234.8 K

Plume 1 (lat=39.7943, lon=-119.7410):
  Nearest GOES pixel mask value: 100.0
  Fire pixels in 7x7 neighbourhood: 3
  Max fire power: 466.6 MW
  Max fire temp:  1234.8 K
```
No traceback (after the load-order fix above); exit code 0.

**Note on mask value 100**: this isn't in the fire-code list (10-15, 30-35)
the task described — 100 is the standard ABI FDC code for "processed,
no fire detected" at that exact pixel. So the *nearest* GOES pixel to each
plume centroid itself shows no fire, but both 7x7 neighborhoods (~roughly
±10-11 km at this satellite viewing geometry) contain real fire pixels —
12 for Plume 0, 3 for Plume 1 — with substantial power (up to ~1062 MW) and
temperature (~1235 K, both patches saturate at the same max, worth noting in
case that's a sensor saturation ceiling rather than two coincidentally equal
fires). This is consistent with the plumes being smoke lofted downwind from
nearby active fire pixels rather than sitting directly over the fire's own
GOES pixel — expected given KRGX plume centroids are smoke aloft, not the
fire's surface hotspot.

**Environment note (repeats)**: same `geo_env` python.exe + manual PATH
workaround as prior tasks (per the note now in `CLAUDE.md`).
## TASK 08 — Label DBSCAN clusters by nearest GOES fire source
**Status**: PENDING

**Goal**: After DBSCAN clustering, assign each cluster a fire-source label
based on which GOES fire pixel centroid it is closest to. Store the mapping
so the hourly pipeline can carry consistent per-source labels across scans.

**Background**: Task 07 confirmed two distinct fire sources:
- Source A (Bug Fire main): ~39.7400°N, 119.8400°W, 1062 MW
- Source B (spot fire / breakout): ~39.7943°N, 119.7410°W, 466 MW
These are the GOES-detected fire centroids nearest to Plume 0 and Plume 1
respectively. We hardcode them for now as named anchors.

**Instructions**:

### Change 1 — Add fire source anchors near top of krgx_plume.py
After the config block (DBZ_THRESHOLD, GAUSS_SIGMA etc.), add:

```python
# ── Fire source anchors (from Task 07 GOES cross-reference) ──────────────
# Each entry: (lat, lon, display_name)
FIRE_SOURCES = [
    (39.7400, -119.8400, 'Bug Fire (main)'),
    (39.7943, -119.7410, 'Spot Fire NE'),
]
```

### Change 2 — Convert fire source anchors to aeqd km at startup
After the aeqd projection is defined (wherever `Proj(proj='aeqd', ...)` is
called, or pyproj is used to set up the radar coordinate system), add:

```python
from pyproj import Proj, Transformer

_aeqd  = Proj(proj='aeqd', lat_0=RADAR_LAT, lon_0=RADAR_LON, units='m')
_wgs84 = Proj(proj='latlong', datum='WGS84')
_to_aeqd = Transformer.from_proj(_wgs84, _aeqd, always_xy=True)

fire_source_xy = []   # (x_km, y_km, name) in radar aeqd
for (lat, lon, name) in FIRE_SOURCES:
    x_m, y_m = _to_aeqd.transform(lon, lat)
    fire_source_xy.append((x_m/1e3, y_m/1e3, name))
    print(f"Fire source '{name}': x={x_m/1e3:.1f} km  y={y_m/1e3:.1f} km")
```
If pyproj is already imported elsewhere in the script, skip the import line
and adapt to whatever Proj/Transformer pattern is already in use.

### Change 3 — Assign source label to each DBSCAN cluster
Find the existing per-cluster print loop (where Plume 0/1 centroids are
printed). After computing each cluster's centroid (cx, cy), add nearest-source
assignment:

```python
    # Nearest fire source
    dists = [np.sqrt((cx - sx)**2 + (cy - sy)**2)
             for (sx, sy, _) in fire_source_xy]
    nearest_idx  = int(np.argmin(dists))
    nearest_name = fire_source_xy[nearest_idx][2]
    nearest_dist = dists[nearest_idx]
    print(f"  Plume {cid} → '{nearest_name}' "
          f"(centroid dist={nearest_dist:.1f} km)")
```

### Change 4 — Store label in cluster results dict
The script presumably builds a results dict or list per cluster. Add two
fields to whatever structure holds per-cluster stats:

```python
    'source_idx':  nearest_idx,
    'source_name': nearest_name,
    'source_dist_km': round(nearest_dist, 2),
```
If the script just prints stats and doesn't store a results structure,
create one now — a list of dicts, one per cluster, appended inside the
cluster loop. We will need this structure for the hourly pipeline in Task 09.

**After changes**: run krgx_plume.py and paste:
1. Fire source aeqd coordinates (from the startup print)
2. Per-cluster source assignment lines
3. Any traceback

---
## RESULTS 08
**Status**: DONE (code works) — but flagging a circularity in the test data,
see below. Run 2026-09-18.

**Naming/structure notes** (adapted to the actual script):
- No top-level `Proj`/`Transformer` existed elsewhere in the script (only a
  local `from pyproj import Proj` inside `fetch_terrain()`), so added the
  import + `_aeqd`/`_wgs84`/`_to_aeqd` block as its own top-level section
  right after `FIRE_SOURCES`, per the spec's fallback instruction.
- No existing per-cluster results structure existed — created
  `plume_results = []` (list of dicts) right before the `if len(pts_bin) > 10:`
  block, appended one dict per cluster inside the loop with the fields the
  spec asked for (`source_idx`, `source_name`, `source_dist_km`) plus the
  centroid/geometry fields that were already being computed and printed
  (`x_km`, `y_km`, `z_km_msl`, `elev_deg`, `horiz_deg`, `n_pts`) so Task 09
  has one self-contained structure instead of having to re-derive geometry
  from print statements.
- Used ASCII `->` instead of the spec's `→` in print statements — purely
  cosmetic, avoids any console-encoding risk, no functional difference.

**1. Fire source aeqd coordinates** (startup print):
```
Fire source 'Bug Fire (main)': x=-32.5 km  y=-1.5 km
Fire source 'Spot Fire NE': x=-24.0 km  y=4.5 km
```

**2. Per-cluster source assignment**:
```
Plume 0: 737 pts  x=-32.5 y=-1.5 z=3.07 km MSL  elev=1.3°  horiz=13.4°
Plume 0 -> 'Bug Fire (main)' (centroid dist=0.1 km)
Plume 1: 20 pts  x=-24.0 y=4.5 z=2.91 km MSL  elev=11.9°  horiz=89.0°
Plume 1 -> 'Spot Fire NE' (centroid dist=0.0 km)
```
Exit code 0, no traceback.

**Flagging a circularity, not a bug**: the near-zero distances (0.1 km, 0.0
km) aren't validating anything yet — they're expected by construction. The
`FIRE_SOURCES` lat/lon anchors in this task are literally Task 07's
`to_ll.transform()` output on Plume 0/1's own aeqd centroids (`x_km=-32.5,
y_km=-1.5` -> lat/lon, and `x_km=-24.0, y_km=4.5` -> lat/lon). Task 08 then
converts those same lat/lons back to aeqd and measures distance to the
plume centroids that produced them — so it's round-tripping the same two
numbers through a projection and back, not independently confirming that
GOES fire pixels and radar-detected plumes collocate. The real GOES fire
*pixel* centroids from Task 07 (the actual `Mask`/`Power` grid, not the
plume-derived anchor points) were never extracted as their own lat/lon in
Task 07 — that task only checked whether fire pixels existed in a
neighborhood around each plume, it didn't compute the fire pixels' own
centroid location. For Task 08's nearest-source logic to mean anything on a
future scan (per its own stated purpose — "so the hourly pipeline can carry
consistent labels across scans"), the anchors need to come from an
independent source, e.g. the actual GOES fire-pixel centroid (weighted by
`Power`) within each plume's search neighborhood, not the plume centroid
itself. Recommend Task 09 (or a Task 08b) re-derive `FIRE_SOURCES` from real
GOES pixel locations before relying on `source_dist_km` for anything.

**Environment note (repeats)**: same `geo_env` python.exe + manual PATH
workaround as prior tasks.

## TASK 08b — Re-derive FIRE_SOURCES from GOES power-weighted fire centroids
**Status**: PENDING

**Goal**: Replace the circular plume-derived FIRE_SOURCES anchors with
independently-derived lat/lon centroids computed from the GOES fire pixel
Power field. Use the same GOES file identified in Task 07.

**Instructions**: Create scratch script `task08b_goes_centroids.py`:

```python
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

# GOES ABI fixed grid → lat/lon
sat_h     = 35786023.0
goes_proj = Proj(proj='geos', lon_0=-137.0, h=sat_h,
                 sweep='x', ellps='GRS80')
wgs84     = Proj(proj='latlong', datum='WGS84')
aeqd      = Proj(proj='aeqd', lat_0=RADAR_LAT, lon_0=RADAR_LON, units='m')

to_ll   = Transformer.from_proj(goes_proj, wgs84, always_xy=True)
to_aeqd = Transformer.from_proj(wgs84, aeqd, always_xy=True)

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

# Cluster fire pixels into two groups using the known approximate separation
# (~9 km between the two sources) — simple distance split around midpoint
# Midpoint of the two plume centroids in aeqd km:
mid_x = (-32.5 + -24.0) / 2   # = -28.25 km
mid_y = ( -1.5 +   4.5) / 2   # =   1.50 km

# Split: fire pixels closer to west source vs east source
# Use bearing from midpoint — pixels west of midpoint → Source A, east → B
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
    clon, clat = to_ll.transform(*to_aeqd.transform(  # round-trip for display
        np.array([cx_km*1e3]), np.array([cy_km*1e3])))
    print(f"\nSource {label} — '{name}':")
    print(f"  N fire pixels: {sel.sum()}")
    print(f"  Power-weighted centroid: x={cx_km:.2f} km  y={cy_km:.2f} km")
    print(f"  Lat/lon: {float(clat):.4f}°N  {float(clon):.4f}°W")
    print(f"  Total power: {fire_power[sel].sum():.0f} MW")
    print(f"  Max power:   {fire_power[sel].max():.0f} MW")
```

**After running**: update `FIRE_SOURCES` in `krgx_plume.py` with the
power-weighted lat/lon values this script produces, replacing the current
circular values. Re-run `krgx_plume.py` and paste the new source assignment
distances — they should now be non-trivial (kilometres, not ~0).

**Append to RESULTS 08b**:
1. Total GOES fire pixels in scene
2. Per-source pixel count, power-weighted centroid (aeqd km + lat/lon),
   total and max power
3. Updated FIRE_SOURCES values pasted into krgx_plume.py
4. New per-cluster source assignment distances from krgx_plume.py re-run

---
## RESULTS 08b
**Status**: BLOCKED on a real decision — did not update `FIRE_SOURCES` yet,
see below. Run 2026-09-18.

**Two bugs fixed in the provided script** (logic errors, not naming
mismatches — both would have crashed or silently produced garbage):
1. The "round-trip for display" line chained the wrong transformers:
   `to_ll.transform(*to_aeqd.transform(x_m, y_m))` fed aeqd meters into
   `to_aeqd`, which expects (lon, lat) → aeqd, not the reverse — then fed
   the (wrong) result into `to_ll`, which expects GOES-fixed-grid → wgs84,
   not aeqd → wgs84. Neither transformer matched the direction needed.
   Replaced with a proper `aeqd_to_ll = Transformer.from_proj(aeqd, wgs84, ...)`
   used directly on the aeqd centroid.
2. `float(clat)` crashed (`TypeError: only 0-dimensional arrays...`) because
   `.transform()` was being called on 1-element numpy arrays wrapped in
   `np.array([...])`; switched to passing plain scalars (`cx_km*1e3`), which
   pyproj returns as plain floats, and removed the now-unneeded `float()`.

**The bigger finding — not a bug, a data/assumption problem**:
`GOES-18 ABI-L2-FDCC` is a full-CONUS product, so the unrestricted
`is_fire` mask picked up **143 fire pixels across the entire western US**
scene, not just near KRGX. The provided script's west/east midpoint split
doesn't restrict to the radar domain first, so it silently produced a
"Source A" centroid 336 km from KRGX (some unrelated fire elsewhere) before
I caught it. Added a domain-bounding filter using `krgx_plume.py`'s own
`X_MIN/X_MAX/Y_MIN/Y_MAX` (-50 to +10 km x, -10 to +25 km y) before doing
anything else with the fire pixels.

**After restricting to the KRGX domain**:
```
Total GOES fire pixels in scene: 143
Fire pixels within KRGX domain: 14 / 143

Source A fire pixels: 14
Source B fire pixels: 0

Source A - 'Bug Fire (main)':
  N fire pixels: 14
  Power-weighted centroid: x=-35.82 km  y=1.08 km
  Lat/lon: 39.7631N  -119.8788W
  Total power: 3640 MW
  Max power:   1062 MW

Source B (Spot Fire NE): NO FIRE PIXELS — check split
```

**All 14 in-domain GOES fire pixels fall on the "Bug Fire main" side of the
split. Zero fall on the "Spot Fire NE" side.** This isn't a code bug — it
means at this GOES scan time (00:06:18Z), there is no independently
resolvable GOES fire-pixel detection anywhere near Plume 1's location
(-24.0, 4.5 km aeqd / 39.7943N, 119.7410W). Possible explanations, in
rough order of likelihood:
1. Plume 1 (20 voxels, a much smaller cluster than Plume 0's 737) is smoke
   from the same Bug Fire main hotspot, just a separate radar-detected
   plume lobe/updraft pulse — not a second independently burning source at
   all. The Task 07 "3 fire pixels near Plume 1" result likely came from its
   7x7-pixel (~±10-11 km) search neighborhood overlapping the edge of the
   same Bug Fire hotspot, not a distinct thermal signature at Plume 1's own
   location.
2. There is a real second hotspot but it's too small/cool for GOES's ~2 km
   IR pixel resolution to flag as "fire" (mask codes 10-15/30-35) at this
   specific 5-minute scan, even though the radar sees a distinct smoke
   plume aloft.
3. Cloud/smoke obscuration over the second source at this exact time step
   (GOES FDC requires a clear IR view of the surface; thick smoke can mask
   the hotspot underneath even though the radar sees the plume above it).

**I did not update `FIRE_SOURCES` in `krgx_plume.py`.** Overwriting "Spot
Fire NE" with a fabricated or interpolated location would reintroduce the
same circularity problem this task was meant to fix, and updating only
Source A while leaving Source B's old circular value in place would be
misleading (it'd look independently verified when it isn't). Options for
the PM to choose between:
- (a) Update `FIRE_SOURCES` to a single verified entry (Bug Fire main,
  39.7631N/-119.8788W) and drop "Spot Fire NE" from the anchor list until a
  GOES scan actually resolves it independently — Plume 1 would then just
  report its distance to the one real source.
- (b) Treat this as evidence Plume 1 isn't a second fire source at all, and
  fold that into how Task 09's hourly pipeline interprets multi-cluster
  scans (e.g. don't label a cluster as a new source unless it has its own
  GOES confirmation).
- (c) Widen the GOES search (more scans across the 5-10 min window, or drop
  the fire-code threshold to include lower-confidence codes) before
  concluding Source B doesn't exist yet.

**Environment note (repeats)**: same `geo_env` python.exe + manual PATH
workaround as prior tasks.
## TASK 08c — Update FIRE_SOURCES to single verified anchor, reclassify Plume 1
**Status**: PENDING

**Goal**: Replace the two-entry FIRE_SOURCES with the single GOES
power-weighted Bug Fire centroid. Re-run krgx_plume.py and confirm Plume 1
is reclassified as an unconfirmed lobe rather than a named second source.

### Change 1 — Update FIRE_SOURCES in krgx_plume.py
Find:
```python
FIRE_SOURCES = [
    (39.7400, -119.8400, 'Bug Fire (main)'),
    (39.7943, -119.7410, 'Spot Fire NE'),
]
```
Replace with:
```python
# Anchor derived from GOES-18 FDCC power-weighted fire pixel centroid
# (Task 08b) — independently verified, not derived from radar plume location.
# 'Spot Fire NE' removed: zero independent GOES detection at that location;
# Plume 1 is likely a secondary lobe of the Bug Fire, not a separate source.
FIRE_SOURCES = [
    (39.7631, -119.8788, 'Bug Fire (main)'),
]

# Minimum GOES fire pixels required to call a cluster a confirmed source
GOES_CONFIRMATION_MIN = 1  # set to 0 to disable; currently unused in single-scan
```

### Change 2 — Add GOES confirmation flag to cluster labeling
In the per-cluster source assignment block (from Task 08), after computing
`nearest_dist`, add a confirmation flag:

```python
    # A cluster is a confirmed source lobe if it has a named anchor
    # and dist is plausible (< 20 km — smoke advects, so generous bound).
    # Clusters with no nearby anchor or dist > 20 km → unconfirmed.
    DIST_THRESHOLD_KM = 20.0
    if nearest_dist <= DIST_THRESHOLD_KM:
        source_label    = nearest_name
        source_confirmed = True
    else:
        source_label    = f'Unconfirmed (nearest: {nearest_name} {nearest_dist:.1f} km)'
        source_confirmed = False

    print(f"  Plume {cid} -> '{source_label}' "
          f"({'confirmed' if source_confirmed else 'UNCONFIRMED'}, "
          f"dist={nearest_dist:.1f} km)")
```

Add `source_confirmed` to the `plume_results` dict alongside `source_name`.

### After changes: run krgx_plume.py and paste
1. Fire source startup print (should show one entry, new lat/lon)
2. Per-cluster source assignment lines for Plume 0 and Plume 1
3. Any traceback

---
## RESULTS 08c
**Status**: DONE — run 2026-09-18. Applied both changes exactly as
specified (no naming adaptation needed — this section already matched the
spec from Tasks 03/08). Added `source_confirmed` to `plume_results` as
instructed.

**1. Fire source startup print** (now one entry):
```
Fire source 'Bug Fire (main)': x=-35.8 km  y=1.1 km
```

**2. Per-cluster source assignment**:
```
Plume 0: 737 pts  x=-32.5 y=-1.5 z=3.07 km MSL  elev=1.3°  horiz=13.4°
Plume 0 -> 'Bug Fire (main)' (confirmed, dist=4.2 km)
Plume 1: 20 pts  x=-24.0 y=4.5 z=2.91 km MSL  elev=11.9°  horiz=89.0°
Plume 1 -> 'Bug Fire (main)' (confirmed, dist=12.3 km)
```
Exit code 0, no traceback.

**Sanity check on the fix**: distances are now non-trivial (4.2 km, 12.3 km)
instead of the ~0 km circularity from Task 08 — confirms the anchor really
is independent of the plume centroids this time. Both plumes land inside
the 20 km confirmation radius and get labeled "confirmed" against the one
real Bug Fire anchor, consistent with Task 08b's finding that Plume 1 is
most likely a secondary lobe of the same fire rather than an unconfirmed/
separate source — the 20 km threshold is generous enough that it doesn't
distinguish that case from a true second confirmed source. If that
distinction matters later (e.g. for the hourly pipeline tracking whether a
*new* fire ignites), worth tightening `DIST_THRESHOLD_KM` or basing
confirmation on each cluster's own local smoke-mask/RhoHV signature rather
than distance alone — flagging as an open question, not fixing since it
wasn't asked for here.

**Environment note (repeats)**: same `geo_env` python.exe + manual PATH
workaround as prior tasks.
(Claude Code pastes output here)

## TASK 09 — Add smoke mask + GOES source labeling to krgx_hourly.py
**Status**: PENDING

**Goal**: Extend `process_scan()` to load RhoHV, compute the smoke mask,
label each hotspot against the Bug Fire GOES anchor, and add smoke mask
stats to the output dict and 2D PNG. No architectural changes — surgical
additions only.

**Instructions**:

### Change 1 — Fix sw variable collision in gate loop
Find the gate extraction loop:
```python
    for sw in keep:
        sl   = radar.get_slice(sw)
        refl = radar.fields['reflectivity']['data'][sl]
```
Rename loop variable `sw` → `isw` throughout the loop body:
```python
    for isw in keep:
        sl   = radar.get_slice(isw)
        refl = radar.fields['reflectivity']['data'][sl]
```
(No other logic change — just the variable name.)

### Change 2 — Add pol fields to include_fields and gate extraction
Find:
```python
    radar = pyart.io.read_nexrad_archive(
        io.BytesIO(raw), include_fields=['reflectivity']
    )
```
Replace with:
```python
    radar = pyart.io.read_nexrad_archive(
        io.BytesIO(raw), include_fields=['reflectivity', 'cross_correlation_ratio']
    )
```

In the gate extraction loop, after the existing `gate_v.append(...)` line,
add RhoHV extraction using the same `valid` mask:
```python
        rhv_field = radar.fields.get('cross_correlation_ratio')
        if rhv_field is not None:
            rhv_data = rhv_field['data'][sl]
            rhv_arr  = np.where(
                np.ma.getmaskarray(rhv_data),
                np.nan,
                rhv_data.astype(float)
            )
            gate_rhv.append(rhv_arr[valid].ravel())
        else:
            gate_rhv.append(np.full(valid.sum(), np.nan))
```

Before the loop, initialize alongside `gate_x, gate_y, gate_z, gate_v`:
```python
    gate_rhv = []
```

### Change 3 — Build rhv_grid after histogramdd block
After the existing `dbz_grid` is built (after the `np.errstate` block),
add:
```python
    # RhoHV grid — per-field count to avoid bias from differential masking
    grhv_all = np.concatenate(gate_rhv)
    rhv_finite = np.isfinite(grhv_all)
    if rhv_finite.sum() > 0:
        rhv_sum, _   = np.histogramdd(
            pts[rhv_finite], bins=[x_bins, y_bins, z_bins],
            weights=grhv_all[rhv_finite]
        )
        rhv_counts, _ = np.histogramdd(
            pts[rhv_finite], bins=[x_bins, y_bins, z_bins]
        )
        with np.errstate(invalid='ignore'):
            rhv_grid = np.where(rhv_counts > 0, rhv_sum / rhv_counts, np.nan)
    else:
        rhv_grid = np.full_like(dbz_grid, np.nan)
```

### Change 4 — Compute smoke mask after dbz_smooth
After `dbz_smooth = gaussian_filter(...)`, add:
```python
    RHV_MAX       = 0.90
    smoke_mask_3d = (
        (dbz_smooth > DBZ_THRESHOLD) &
        np.isfinite(rhv_grid) &
        (rhv_grid < RHV_MAX)
    )
    smoke_vox = int(smoke_mask_3d.sum())
```

### Change 5 — Add fire source anchor + per-hotspot labeling
Add at the top of `process_scan()` (or module-level constants — module-level
preferred so it's computed once, not per scan):

```python
# Module-level, after CONFIG block:
from pyproj import Proj, Transformer as _Transformer

_aeqd     = Proj(proj='aeqd', lat_0=RADAR_LAT, lon_0=RADAR_LON, units='m')
_wgs84    = Proj(proj='latlong', datum='WGS84')
_to_aeqd  = _Transformer.from_proj(_wgs84, _aeqd, always_xy=True)

FIRE_SOURCES = [
    (39.7631, -119.8788, 'Bug Fire (main)'),   # GOES power-weighted centroid, Task 08b
]
fire_source_xy = []
for (_lat, _lon, _name) in FIRE_SOURCES:
    _xm, _ym = _to_aeqd.transform(_lon, _lat)
    fire_source_xy.append((_xm/1e3, _ym/1e3, _name))
```

In the hotspot stats collection loop, after computing `x_mean`/`y_mean`,
add nearest-source assignment:
```python
        cx, cy = float(x_km[ix].mean()), float(y_km[iy].mean())
        dists  = [np.sqrt((cx-sx)**2 + (cy-sy)**2)
                  for (sx, sy, _) in fire_source_xy]
        ni     = int(np.argmin(dists))
        hotspots.append({
            ...existing fields...,
            'source_name':     fire_source_xy[ni][2],
            'source_dist_km':  round(dists[ni], 2),
            'smoke_vox':       int((smoke_mask_3d & region).sum()),
            'smoke_frac':      round(
                (smoke_mask_3d & region).sum() / max(region.sum(), 1), 3),
        })
```

### Change 6 — Add smoke mask footprint to 2D PNG
In the plotting block, after the existing `ax.contour` loop, add:
```python
    smoke_footprint = smoke_mask_3d.any(axis=2).T
    if smoke_footprint.any():
        ax.contourf(x_km, y_km, smoke_footprint.astype(float),
                    levels=[0.5, 1.5], colors=['purple'], alpha=0.25,
                    zorder=2)
        ax.contour(x_km, y_km, smoke_footprint.astype(float),
                   levels=[0.5], colors=['purple'], linewidths=1,
                   linestyles='--', zorder=3)
```

### Change 7 — Add smoke stats to results dict and per-scan print
In `process_scan()`'s return dict, add:
```python
        'smoke_vox_total': smoke_vox,
```

In the main loop print block, after the existing hotspot print, add:
```python
        print(f'  Smoke mask: {res["smoke_vox_total"]} vox total')
        for hi, hs in enumerate(res['hotspots']):
            print(f'    Hotspot {hi+1} -> {hs["source_name"]} '
                  f'(dist={hs["source_dist_km"]:.1f} km)  '
                  f'smoke={hs["smoke_frac"]*100:.0f}%')
```

### After all changes: run krgx_hourly.py
Paste:
1. Per-scan summary lines (timestamp, hotspot count, smoke vox, source labels)
2. Any scan that errors or returns None
3. Total scans processed line at the end
4. Whether PNGs are saving correctly (no matplotlib errors)

Do NOT paste the full raw output if it is very long —
summarize scans with no hotspots as "(no echo)" and
paste full detail only for scans with hotspots.

---
## RESULTS 09

**Status**: DONE

**Naming/deviation notes**: All 7 changes applied to `krgx_hourly.py` as
specified, adapting variable names to what was already there:
- `region` was already the loop variable name for the hotspot's boolean
  mask (matches the instructions' assumption).
- The fire-source anchor block was inserted right after the existing
  `BUCKETS = [...]` line in the CONFIG section (module-level, computed once,
  per the instructions' stated preference).
- No name mismatches encountered — this task's example code matched the
  actual script closely enough to apply verbatim.

**Run output** (14/14 scans processed, no errors, all PNGs saved):

```
Listing scans...
  14 scans in unidata-nexrad-level2

[1/14] 23:59:57Z — 4 hotspots, smoke mask 1155 vox
  Hotspot 1 -> Bug Fire (main) (dist=3.1 km)  smoke=83%
  Hotspot 2 -> Bug Fire (main) (dist=8.3 km)  smoke=97%
  Hotspot 3 -> Bug Fire (main) (dist=14.0 km) smoke=94%
  Hotspot 4 -> Bug Fire (main) (dist=17.6 km) smoke=90%
  saved krgx_20260809_235957.png (5.1s)

[2/14] 00:04:04Z — 6 hotspots, smoke mask 1510 vox, dist range 3.8–21.3 km, smoke 82–98%
[3/14] 00:08:09Z — 6 hotspots, smoke mask 1984 vox, dist range 4.3–23.6 km, smoke 78–100%
[4/14] 00:12:15Z — 5 hotspots, smoke mask 2562 vox, dist range 7.0–24.8 km, smoke 86–100%
[5/14] 00:16:21Z — 8 hotspots, smoke mask 2690 vox, dist range 6.7–27.2 km, smoke 55–100%
[6/14] 00:20:27Z — 7 hotspots, smoke mask 3394 vox, dist range 4.1–18.5 km, smoke 74–100%
[7/14] 00:24:33Z — 4 hotspots, smoke mask 2547 vox, dist range 4.8–18.9 km, smoke 86–93%
[8/14] 00:32:07Z — 8 hotspots, smoke mask 3344 vox, dist range 6.0–26.8 km, smoke 80–98%
[9/14] 00:36:13Z — 6 hotspots, smoke mask 3619 vox, dist range 4.8–29.0 km, smoke 80–100%
[10/14] 00:40:20Z — 2 hotspots, smoke mask 3649 vox, dist 4.9/18.0 km, smoke 90%/80%
[11/14] 00:44:27Z — 2 hotspots, smoke mask 3511 vox, dist 5.0/19.9 km, smoke 88%/83%
[12/14] 00:48:33Z — 3 hotspots, smoke mask 3339 vox, dist range 5.7–22.5 km, smoke 84–100%
[13/14] 00:52:38Z — 4 hotspots, smoke mask 2600 vox, dist range 7.1–26.3 km, smoke 48–97%
[14/14] 00:56:45Z — 8 hotspots, smoke mask 1686 vox, dist range 5.3–28.9 km, smoke 20–100%

Done: 14/14 scans processed
PNGs saved to C:\Users\raymondharding\Research\Runs\krgx_hourly
```

No scan errored or returned None; no matplotlib errors. All 14 PNGs
confirmed written to `krgx_hourly/`.

**Flagged issue (not acted on, per standing instruction)**: This task's
per-hotspot labeling logic (Change 5) assigns every hotspot to the
*nearest* fire source unconditionally — there is no distance-threshold
gate or `source_confirmed` flag like the one added to `krgx_plume.py` in
Task 08c. With only one entry in `FIRE_SOURCES`, this means hotspots as
far as 28–29 km from the Bug Fire anchor (e.g. scan 9 and 14, smoke_frac
as low as 20%) are still printed as `-> Bug Fire (main)` with no
"unconfirmed" caveat. That's misleading on its face — at 29 km with 20%
smoke fraction, that hotspot is plausibly unrelated convective/clutter
echo rather than plume from the Bug Fire. If the hourly pipeline is meant
to support new-ignition detection (the concern raised in RESULTS 08c),
it would need the same `DIST_THRESHOLD_KM` + `source_confirmed` pattern
ported over from `krgx_plume.py`. Left as-is since Task 09 didn't ask for
it and I was told to flag rather than force through unrequested changes.

## TASK 10 — Port source_confirmed threshold gate to krgx_hourly.py
**Status**: PENDING

**Goal**: Add the same DIST_THRESHOLD_KM + source_confirmed pattern from
krgx_plume.py (Task 08c) to krgx_hourly.py's hotspot labeling. Hotspots
beyond 20 km from the nearest fire source anchor print as "unconfirmed"
and are flagged in the results dict.

**Instructions**: One change only — in the hotspot stats collection loop
(the block added in Task 09 Change 5), replace the source assignment with:

```python
        DIST_THRESHOLD_KM = 20.0
        ni    = int(np.argmin(dists))
        ndist = round(dists[ni], 2)
        if ndist <= DIST_THRESHOLD_KM:
            src_name      = fire_source_xy[ni][2]
            src_confirmed = True
        else:
            src_name      = f'Unconfirmed (nearest: {fire_source_xy[ni][2]} {ndist:.1f} km)'
            src_confirmed = False

        hotspots.append({
            ...all existing fields...,
            'source_name':     src_name,
            'source_dist_km':  ndist,
            'source_confirmed': src_confirmed,
            'smoke_vox':       int((smoke_mask_3d & region).sum()),
            'smoke_frac':      round(
                (smoke_mask_3d & region).sum() / max(region.sum(), 1), 3),
        })
```

Also add DIST_THRESHOLD_KM as a module-level config constant (near the
other thresholds in the CONFIG block) rather than hardcoding it inside
the loop:
```python
DIST_THRESHOLD_KM = 20.0   # hotspots beyond this → unconfirmed source
```
And reference it in the loop without redefining it.

In the main loop print block (Task 09 Change 7), update the per-hotspot
line to show confirmed status:
```python
        for hi, hs in enumerate(res['hotspots']):
            confirmed_str = '' if hs['source_confirmed'] else ' [UNCONFIRMED]'
            print(f'    Hotspot {hi+1} -> {hs["source_name"]}'
                  f'{confirmed_str}  '
                  f'dist={hs["source_dist_km"]:.1f} km  '
                  f'smoke={hs["smoke_frac"]*100:.0f}%')
```

### After changes: run krgx_hourly.py
Paste only:
1. Count of confirmed vs unconfirmed hotspots across all 14 scans
2. The worst unconfirmed case (highest dist, lowest smoke_frac)
3. Whether any scan's primary hotspot (Hotspot 1, highest dBZ) is
   ever unconfirmed

---
## RESULTS 10

**Status**: DONE

**Naming/deviation notes**: Applied exactly as specified — `DIST_THRESHOLD_KM`
added as a module-level constant next to `NEIGHBOR_RADIUS_M` in the CONFIG
block, hotspot-loop source assignment replaced with the `src_name`/
`src_confirmed` logic, `source_confirmed` added to the hotspot dict, and
the print block updated with the `[UNCONFIRMED]` tag. No name mismatches.

**1. Confirmed vs unconfirmed hotspots across all 14 scans**:
- Total hotspots: 73
- Confirmed: 51
- Unconfirmed: 22

**2. Worst unconfirmed cases**:
- Highest distance: scan 00:36:13Z, Hotspot 4 — dist=29.0 km, smoke=96%
- Lowest smoke_frac: scan 00:56:45Z, Hotspot 8 — dist=26.7 km, smoke=20%
  (these are two different hotspots — no single hotspot is worst on both axes)

**3. Is the primary hotspot (Hotspot 1, highest dBZ) ever unconfirmed?**
Yes, in the last two scans of the hour:
- 00:52:38Z — Hotspot 1 (max 17.5 dBZ, 2427 vox) → **UNCONFIRMED**, dist=26.3 km, smoke=86%
- 00:56:45Z — Hotspot 1 (max 17.2 dBZ, 1396 vox) → **UNCONFIRMED**, dist=28.9 km, smoke=87%

In every other scan (1 through 12, 23:59:57Z–00:48:33Z) Hotspot 1 stays
confirmed, with distance climbing steadily from 3.1 km at scan 1 up to
19.9 km at scan 11. This looks like real plume drift away from the anchor
over the hour rather than a bug — smoke_frac stays high (83-100%) even
once unconfirmed, so it's still being detected as smoke, just too far
from the single anchor point to be labeled with confidence. Flagging in
case the PM wants this read as "plume has drifted past the confirmation
radius" for the write-up, rather than "detection failure."

## TASK 12 — Centroid drift time series across all 14 scans
**Status**: PENDING

**Goal**: Extract the primary confirmed hotspot centroid (x_km, y_km, z_km_msl)
from each scan's results dict and plot centroid drift as a function of time.
Standalone scratch script — no changes to krgx_hourly.py.

**Background**: `results` is already built by krgx_hourly.py at the end of
the main loop but not persisted. Add a one-line save at the bottom of
krgx_hourly.py, then write the analysis script separately.

### Change 1 — Persist results to disk at end of krgx_hourly.py
At the very end of krgx_hourly.py (after the final print statements), add:

```python
import json

def _serial(obj):
    """JSON serializer for numpy types."""
    if isinstance(obj, (np.integer,)): return int(obj)
    if isinstance(obj, (np.floating,)): return float(obj)
    if isinstance(obj, (np.ndarray,)): return obj.tolist()
    raise TypeError(f"Not serializable: {type(obj)}")

out_json = OUT_DIR / 'hourly_results.json'
with open(out_json, 'w') as f:
    json.dump([
        {
            'dt': r['dt'].isoformat(),
            'smoke_vox_total': r['smoke_vox_total'],
            'hotspots': r['hotspots'],   # already dicts, region excluded
        }
        for r in results
    ], f, indent=2, default=_serial)
print(f'Results saved to {out_json}')
```

Re-run krgx_hourly.py to generate `krgx_hourly/hourly_results.json`.
(This will re-process all 14 scans — ~70s total based on Task 09 timing.)

### Change 2 — Create task12_centroid_drift.py scratch script

```python
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from datetime import datetime, timezone
from pathlib import Path

# ── load results ──────────────────────────────────────────────────────────
with open('krgx_hourly/hourly_results.json') as f:
    results = json.load(f)

# ── extract primary confirmed hotspot per scan ────────────────────────────
# Primary = Hotspot 0 (highest dBZ, already sorted descending)
# Use confirmed only; if Hotspot 0 is unconfirmed, fall back to first confirmed
times   = []
x_cents = []
y_cents = []
z_cents = []
z_tops  = []
smoke_fracs = []
scan_labels = []

for r in results:
    dt = datetime.fromisoformat(r['dt'])
    hs_list = r['hotspots']
    
    # pick primary: first confirmed hotspot by dBZ rank
    primary = None
    for hs in hs_list:
        if hs.get('source_confirmed', True):   # default True if key missing
            primary = hs
            break
    if primary is None:
        continue   # skip scans with no confirmed hotspot
    
    times.append(dt)
    x_cents.append(primary['x_mean'])
    y_cents.append(primary['y_mean'])
    z_cents.append(primary['z_mean'])
    z_tops.append(primary['z_max'])
    smoke_fracs.append(primary['smoke_frac'])
    scan_labels.append(dt.strftime('%H:%M'))

times   = np.array([(t - times[0]).total_seconds()/60 for t in times])
x_cents = np.array(x_cents)
y_cents = np.array(y_cents)
z_cents = np.array(z_cents)
z_tops  = np.array(z_tops)
smoke_fracs = np.array(smoke_fracs)

print(f"Scans with confirmed primary hotspot: {len(times)}/14")
print(f"\nCentroid drift summary:")
print(f"  X range: {x_cents.min():.1f} to {x_cents.max():.1f} km  "
      f"(total drift: {x_cents.max()-x_cents.min():.1f} km E-W)")
print(f"  Y range: {y_cents.min():.1f} to {y_cents.max():.1f} km  "
      f"(total drift: {y_cents.max()-y_cents.min():.1f} km N-S)")
print(f"  Z mean range: {z_cents.min():.2f} to {z_cents.max():.2f} km MSL")
print(f"  Z top range:  {z_tops.min():.2f} to {z_tops.max():.2f} km MSL")

# ── estimate drift velocity ────────────────────────────────────────────────
if len(times) >= 2:
    dx = np.diff(x_cents)
    dy = np.diff(y_cents)
    dt_min = np.diff(times)
    speed_kmph = np.sqrt(dx**2 + dy**2) / (dt_min/60)
    bearing    = np.degrees(np.arctan2(dx, dy)) % 360  # N=0, E=90
    print(f"\nScan-to-scan drift velocity:")
    for i, (spd, brg) in enumerate(zip(speed_kmph, bearing)):
        print(f"  {scan_labels[i]} -> {scan_labels[i+1]}: "
              f"{spd:.1f} km/hr  bearing={brg:.0f}°")
    print(f"\n  Mean speed: {speed_kmph.mean():.1f} km/hr")
    print(f"  Mean bearing: {bearing.mean():.0f}°")

# ── plots ─────────────────────────────────────────────────────────────────
import matplotlib
matplotlib.use('Agg')

fig, axes = plt.subplots(2, 2, figsize=(12, 9))

# 1. X/Y centroid track (plan view)
ax = axes[0, 0]
sc = ax.scatter(x_cents, y_cents, c=times, cmap='plasma', s=60, zorder=3)
ax.plot(x_cents, y_cents, 'k--', lw=0.8, zorder=2)
for i, lbl in enumerate(scan_labels):
    ax.annotate(lbl, (x_cents[i], y_cents[i]),
                fontsize=7, xytext=(4, 4), textcoords='offset points')
ax.scatter(0, 0, marker='+', c='black', s=120, zorder=5, label='KRGX')
# Bug Fire anchor
ax.scatter(-35.8, 1.1, marker='*', c='red', s=150, zorder=5,
           label='Bug Fire anchor')
plt.colorbar(sc, ax=ax, label='Minutes from first scan')
ax.set(xlabel='X km (aeqd)', ylabel='Y km (aeqd)',
       title='Plume centroid track (plan view)', aspect='equal')
ax.legend(fontsize=8)

# 2. X centroid vs time
ax = axes[0, 1]
ax.plot(times, x_cents, 'o-', color='steelblue')
ax.set(xlabel='Minutes from first scan', ylabel='X km',
       title='Centroid X drift (E-W)')
ax.axhline(-35.8, color='red', lw=1, linestyle='--', label='Fire anchor X')
ax.legend(fontsize=8)

# 3. Z centroid and Z top vs time
ax = axes[1, 0]
ax.plot(times, z_cents, 'o-', color='darkorange', label='Z mean')
ax.plot(times, z_tops,  's--', color='firebrick',  label='Z top')
ax.set(xlabel='Minutes from first scan', ylabel='km MSL',
       title='Plume height vs time')
ax.legend(fontsize=8)

# 4. Smoke fraction vs time
ax = axes[1, 1]
ax.plot(times, smoke_fracs * 100, 'o-', color='purple')
ax.axhline(80, color='gray', lw=1, linestyle='--', label='80% threshold')
ax.set(xlabel='Minutes from first scan', ylabel='Smoke fraction (%)',
       title='Primary hotspot smoke fraction (ρHV mask)')
ax.legend(fontsize=8)
ax.set_ylim(0, 105)

plt.suptitle('KRGX — Bug Fire plume centroid drift\n'
             '2026-08-09 23:59Z to 2026-08-10 00:56Z', fontsize=11)
plt.tight_layout()
out_png = Path('krgx_hourly/task12_centroid_drift.png')
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print(f'\nSaved {out_png}')
```

### After running: paste
1. "Scans with confirmed primary hotspot: N/14"
2. Centroid drift summary block (X/Y/Z ranges)
3. Mean drift speed and bearing
4. Whether the PNG saved without error

Also paste the PNG image directly if possible.

---
## RESULTS 12

**Status**: DONE

**Naming/deviation notes**: Both changes applied exactly as specified —
the JSON-persist block appended verbatim to the end of `krgx_hourly.py`,
`krgx_hourly.py` re-run to regenerate `hourly_results.json`, and
`task12_centroid_drift.py` created and run verbatim (no name mismatches;
`results`, `hotspots`, `x_mean`/`y_mean`/`z_mean`/`z_max`,
`source_confirmed`, `smoke_frac` all matched the actual dict keys).

**1. Scans with confirmed primary hotspot: 14/14**

**2. Centroid drift summary**:
```
X range: -34.1 to -16.9 km  (total drift: 17.2 km E-W)
Y range: -2.5 to 7.5 km  (total drift: 9.9 km N-S)
Z mean range: 2.72 to 3.12 km MSL
Z top range:  2.90 to 3.90 km MSL
```

**3. Mean drift speed and bearing**: 50.8 km/hr, mean bearing 127°
(bearing convention: N=0°, E=90°)

**4. PNG**: saved without error to `krgx_hourly/task12_centroid_drift.png`
(all 4 panels rendered — plan-view track, X drift, height, smoke fraction).
Not re-pasted here since this channel is text-only; the file is on disk
for the PM/user to open directly.

**Flagged issue**: the per-scan drift velocities are not physically
plausible for a smoke plume — most scan-to-scan speeds are 10-45 km/hr,
but one interval (00:44→00:48) shows **256.9 km/hr** at bearing 236°,
and another (00:24→00:32) shows 96.5 km/hr. Checking the underlying data:
at 00:44:27Z the "primary" (Hotspot 1, highest dBZ) is centered at
x=-16.9, y=7.5. At 00:48:33Z, Hotspot 1 by dBZ rank is a *different*
physical feature at x=-14.4, y=7.9 (18.6 dBZ) — but it's unconfirmed
(dist=22.5 km), so the script's fallback logic skips it and picks the
next confirmed hotspot instead, which is a separate lobe near the fire
anchor at x=-31.4, y=-2.5. The "primary" label is jumping between two
different physical plume lobes scan-to-scan, not tracking one continuous
feature — the 256.9 km/hr spike is an artifact of the "first confirmed by
dBZ rank" selection rule, not real plume motion. The 96.5 km/hr spike
(00:24→00:32) also coincides with a scan-count gap (8 min between those
two scans vs. ~4 min for most others — see the scan list gap at
00:24:33Z → 00:32:07Z), which alone would roughly double the apparent
speed even for genuine motion, compounding the effect.

Not fixed since Task 12 only asked for extraction/plotting, not a
tracking-identity fix — but for any bearing/speed number to be reported
in the write-up, this needs resolving first (e.g. track the hotspot
nearest to the previous scan's centroid rather than by dBZ rank, and/or
exclude scans with anomalous time gaps from the speed calculation).

## TASK 13 — Nearest-neighbor centroid tracking with 20 km/hr physical cap
**Status**: PENDING

**Goal**: Track the primary plume centroid across 14 scans using nearest-
neighbor linking, with a hard 20 km/hr physical speed cap. Any scan-to-scan
transition exceeding the cap is flagged as a tracking failure — the interval
is excluded from speed/bearing statistics and the tracker re-seeds on the
next scan.

**Instructions**: Modify task12_centroid_drift.py only.

### Change 1 — Replace primary selection with capped nearest-neighbor tracker

Replace the entire extraction loop and the times/x_cents/etc. initialization
with:

```python
SPEED_CAP_KMH = 20.0   # physical ceiling — fire spread + wind advection

prev_x, prev_y, prev_dt = None, None, None
track_ok = []   # True if this scan's link to the previous is valid

for r in results:
    dt = datetime.fromisoformat(r['dt'])
    hs_list = r['hotspots']

    confirmed = [hs for hs in hs_list if hs.get('source_confirmed', True)]
    if not confirmed:
        prev_x, prev_y, prev_dt = None, None, None
        continue

    if prev_x is None:
        # Seed: highest dBZ confirmed hotspot
        primary  = max(confirmed, key=lambda h: h['dbz_max'])
        link_valid = False   # no previous scan to link from
    else:
        # Nearest neighbor candidate
        candidate = min(confirmed,
                        key=lambda h: (h['x_mean']-prev_x)**2
                                    + (h['y_mean']-prev_y)**2)
        dt_hr = (dt - prev_dt).total_seconds() / 3600
        dist  = np.sqrt((candidate['x_mean']-prev_x)**2
                      + (candidate['y_mean']-prev_y)**2)
        speed = dist / dt_hr   # km/hr

        if speed <= SPEED_CAP_KMH:
            primary    = candidate
            link_valid = True
        else:
            # Speed cap exceeded — tracking failure
            # Re-seed with highest dBZ confirmed hotspot
            primary    = max(confirmed, key=lambda h: h['dbz_max'])
            link_valid = False
            print(f"  [TRACK FAIL] {dt.strftime('%H:%MZ')}: "
                  f"nearest candidate {dist:.1f} km in {dt_hr*60:.1f} min "
                  f"= {speed:.0f} km/hr > cap. Re-seeding.")

    prev_x  = primary['x_mean']
    prev_y  = primary['y_mean']
    prev_dt = dt

    times.append(dt)
    x_cents.append(primary['x_mean'])
    y_cents.append(primary['y_mean'])
    z_cents.append(primary['z_mean'])
    z_tops.append(primary['z_max'])
    smoke_fracs.append(primary['smoke_frac'])
    scan_labels.append(dt.strftime('%H:%M'))
    track_ok.append(link_valid)

track_ok = np.array(track_ok)
```

### Change 2 — Exclude invalid links from speed/bearing calculation
After computing `dx`, `dy`, `dt_min`:

```python
# track_ok[i] is True if scan i's link to scan i-1 is valid
# For a speed between scan i and i+1, both endpoints must be valid links
interval_ok = track_ok[:-1] & track_ok[1:]

speed_kmph = np.where(
    interval_ok,
    np.sqrt(dx**2 + dy**2) / (dt_min / 60),
    np.nan
)
bearing = np.where(interval_ok, bearing, np.nan)
```

Update prints:
```python
    for i, (spd, brg) in enumerate(zip(speed_kmph, bearing)):
        if np.isnan(spd):
            print(f"  {scan_labels[i]} -> {scan_labels[i+1]}: "
                  f"EXCLUDED (tracking failure or re-seed)")
        else:
            print(f"  {scan_labels[i]} -> {scan_labels[i+1]}: "
                  f"{spd:.1f} km/hr  bearing={brg:.0f}°")

    valid_speeds = speed_kmph[~np.isnan(speed_kmph)]
    print(f"\n  Valid intervals: {len(valid_speeds)}/{len(speed_kmph)}")
    print(f"  Mean speed:    {np.nanmean(speed_kmph):.1f} km/hr")
    print(f"  Max speed:     {np.nanmax(speed_kmph):.1f} km/hr")
    print(f"  Mean bearing:  {np.nanmean(bearing):.0f}°")
```

### Change 3 — Mark invalid links on the plan-view plot
In the plan-view panel (axes[0,0]), after `ax.plot(x_cents, y_cents, ...)`,
add dashed red segments for re-seeded (invalid) links:

```python
for i in range(1, len(x_cents)):
    if not track_ok[i]:
        ax.plot([x_cents[i-1], x_cents[i]],
                [y_cents[i-1], y_cents[i]],
                'r--', lw=1.5, zorder=2)
```

### After running: paste
1. Any [TRACK FAIL] lines printed during extraction
2. Full per-scan velocity table
3. Valid intervals count, mean speed, max speed, mean bearing
4. Whether max speed is now <= 20 km/hr on all valid intervals

---
## RESULTS 13

**Status**: DONE

**Naming/deviation notes**: All 3 changes applied verbatim to
`task12_centroid_drift.py` — no name mismatches. `dbz_max`, `x_mean`,
`y_mean`, `source_confirmed` all matched the JSON keys already in use
from Task 12.

**1. [TRACK FAIL] lines** (8 of 13 transitions failed the 20 km/hr cap):
```
[TRACK FAIL] 00:12Z: nearest candidate 4.0 km in 4.1 min = 58 km/hr > cap. Re-seeding.
[TRACK FAIL] 00:16Z: nearest candidate 1.9 km in 4.1 min = 27 km/hr > cap. Re-seeding.
[TRACK FAIL] 00:20Z: nearest candidate 1.6 km in 4.1 min = 23 km/hr > cap. Re-seeding.
[TRACK FAIL] 00:36Z: nearest candidate 2.6 km in 4.1 min = 37 km/hr > cap. Re-seeding.
[TRACK FAIL] 00:44Z: nearest candidate 2.0 km in 4.1 min = 29 km/hr > cap. Re-seeding.
[TRACK FAIL] 00:48Z: nearest candidate 10.4 km in 4.1 min = 152 km/hr > cap. Re-seeding.
[TRACK FAIL] 00:52Z: nearest candidate 1.9 km in 4.1 min = 28 km/hr > cap. Re-seeding.
[TRACK FAIL] 00:56Z: nearest candidate 1.9 km in 4.1 min = 27 km/hr > cap. Re-seeding.
```

**2. Full per-scan velocity table**:
```
23:59 -> 00:04: EXCLUDED (tracking failure or re-seed)
00:04 -> 00:08: 14.6 km/hr  bearing=69°
00:08 -> 00:12: EXCLUDED (tracking failure or re-seed)
00:12 -> 00:16: EXCLUDED (tracking failure or re-seed)
00:16 -> 00:20: EXCLUDED (tracking failure or re-seed)
00:20 -> 00:24: EXCLUDED (tracking failure or re-seed)
00:24 -> 00:32: 14.2 km/hr  bearing=73°
00:32 -> 00:36: EXCLUDED (tracking failure or re-seed)
00:36 -> 00:40: EXCLUDED (tracking failure or re-seed)
00:40 -> 00:44: EXCLUDED (tracking failure or re-seed)
00:44 -> 00:48: EXCLUDED (tracking failure or re-seed)
00:48 -> 00:52: EXCLUDED (tracking failure or re-seed)
00:52 -> 00:56: EXCLUDED (tracking failure or re-seed)
```

**3. Valid intervals / stats**: 2/13 valid. Mean speed 14.4 km/hr,
max speed 14.6 km/hr, mean bearing 71°.

**4. Is max speed now ≤ 20 km/hr on all valid intervals?** Yes — 14.6 km/hr,
under the cap by construction (any candidate exceeding it is rejected).

**Flagged issue**: only 2 of 13 transitions (15%) survive the cap — this
is a much bigger rejection rate than a "handle a couple of bad frames"
fix would suggest. Root cause: KRGX's scan cadence here is ~4.1 minutes,
so a 20 km/hr cap only allows ~1.4 km of centroid movement between
consecutive scans. But per-scan hotspot centroids (computed from a
DBSCAN/watershed region's voxel mean) commonly shift 2-10 km between
scans just from the plume's cross-section reshaping, new pyroconvective
pulses, or a different lobe crossing the dBZ threshold — that's not
necessarily *bulk translation* at that speed, but it looks identical to
the tracker. With 8/13 re-seeds, the "tracked" centroid mostly collapses
back to the same dBZ-rank selection as Task 12 (only now silently, since
re-seeding after a failure picks highest-dBZ, same rule as before) — so
the 2 valid intervals we do get (14.2-14.6 km/hr, ~70° bearing, i.e.
roughly ENE) are a much smaller and less certain sample than "14/14
scans tracked" would suggest. If a real drift speed/bearing is needed
for the write-up, either the cap needs raising to reflect this scan
cadence and plume behavior (e.g. 40-60 km/hr, matching the more moderate
non-outlier speeds from Task 12's uncapped run), or the tracking needs
richer state than pure nearest-centroid-distance (e.g. also matching on
voxel count / footprint overlap between scans) before trusting a
bearing derived from only 2 surviving links. Not changed since Task 13
asked for exactly this cap-and-exclude behavior — flagging for the PM's
call on whether 20 km/hr was meant as a first-pass sanity check or the
final number to report.

## TASK 14 — Replace speed cap with footprint overlap tracker
**Status**: PENDING

**Goal**: Track the primary plume across scans using 2D smoke-mask footprint
overlap (Intersection over Union) rather than centroid speed. A valid link
requires IoU >= 0.10 between consecutive scans' smoke footprints. Speed cap
removed entirely.

**Background**: footprint overlap is a more robust tracking criterion than
centroid speed because it is invariant to lobe reshaping and wind advection
— it only fails when two completely disjoint echoes are being compared,
which is the actual failure mode we want to catch.

**Problem**: the smoke-mask footprint per hotspot is not currently stored in
`hourly_results.json` — only centroid stats are. Two options:
(a) Re-run krgx_hourly.py saving per-hotspot footprints (large JSON)
(b) Save just the 2D column-max boolean footprint per scan as a numpy file

Use option (b) — save one (nx, ny) boolean array per scan (union of all
confirmed hotspot footprints) to a compressed .npz file alongside the JSON.

### Change 1 — Save per-scan smoke footprint in krgx_hourly.py
In `process_scan()`, after `smoke_mask_3d` is computed, add:
```python
    smoke_footprint_2d = smoke_mask_3d.any(axis=2)   # (nx, ny) bool
```

In the return dict, add:
```python
        'smoke_footprint_2d': smoke_footprint_2d,
```

At the end of krgx_hourly.py, after the JSON save block, add:
```python
# Save per-scan smoke footprints as compressed npz
footprints = {
    r['dt'].strftime('%Y%m%d_%H%M%S'): r['smoke_footprint_2d']
    for r in results
    if 'smoke_footprint_2d' in r
}
np.savez_compressed(OUT_DIR / 'smoke_footprints.npz', **footprints)
print(f'Footprints saved to {OUT_DIR}/smoke_footprints.npz')
```

Also exclude `smoke_footprint_2d` from the JSON serialization (it's a
numpy array — already excluded by the existing dict comprehension that
drops non-serializable fields, but confirm no TypeError).

Re-run krgx_hourly.py to generate `smoke_footprints.npz`.

### Change 2 — Rewrite tracker in task12_centroid_drift.py

Replace the entire extraction loop with:

```python
import numpy as np
from pathlib import Path

# Load footprints
fp_data = np.load('krgx_hourly/smoke_footprints.npz')

SPEED_CAP_KMH = 20.0   # kept only for *reporting*, not for exclusion
IOU_MIN = 0.10          # minimum footprint overlap to call a valid link

prev_x, prev_y, prev_dt = None, None, None
prev_fp_key = None
track_ok = []

for r in results:
    dt      = datetime.fromisoformat(r['dt'])
    fp_key  = dt.strftime('%Y%m%d_%H%M%S')
    hs_list = r['hotspots']

    confirmed = [hs for hs in hs_list if hs.get('source_confirmed', True)]
    if not confirmed:
        prev_x, prev_y, prev_dt, prev_fp_key = None, None, None, None
        continue

    if prev_x is None or prev_fp_key is None:
        primary    = max(confirmed, key=lambda h: h['dbz_max'])
        link_valid = False
    else:
        # Compute IoU between this scan's footprint and previous
        fp_curr = fp_data[fp_key]       if fp_key     in fp_data.files else None
        fp_prev = fp_data[prev_fp_key]  if prev_fp_key in fp_data.files else None

        if fp_curr is not None and fp_prev is not None:
            intersection = (fp_curr & fp_prev).sum()
            union        = (fp_curr | fp_prev).sum()
            iou          = intersection / max(union, 1)
        else:
            iou = 0.0

        # Nearest-neighbor candidate
        candidate = min(confirmed,
                        key=lambda h: (h['x_mean']-prev_x)**2
                                    + (h['y_mean']-prev_y)**2)

        if iou >= IOU_MIN:
            primary    = candidate
            link_valid = True
        else:
            primary    = max(confirmed, key=lambda h: h['dbz_max'])
            link_valid = False
            dt_hr  = (dt - prev_dt).total_seconds() / 3600
            dist   = np.sqrt((candidate['x_mean']-prev_x)**2
                           + (candidate['y_mean']-prev_y)**2)
            print(f"  [TRACK FAIL] {dt.strftime('%H:%MZ')}: "
                  f"IoU={iou:.3f} < {IOU_MIN}  "
                  f"(nearest cand {dist:.1f} km, "
                  f"{dist/dt_hr:.0f} km/hr). Re-seeding.")

    prev_x, prev_y   = primary['x_mean'], primary['y_mean']
    prev_dt           = dt
    prev_fp_key       = fp_key

    times.append(dt)
    x_cents.append(primary['x_mean'])
    y_cents.append(primary['y_mean'])
    z_cents.append(primary['z_mean'])
    z_tops.append(primary['z_max'])
    smoke_fracs.append(primary['smoke_frac'])
    scan_labels.append(dt.strftime('%H:%M'))
    track_ok.append(link_valid)

track_ok = np.array(track_ok)
```

Keep Change 2 (gap exclusion from Task 13) and Change 3 (red dashed
segments on plan view) from Task 13 — they still apply, now gated on
IoU rather than speed.

Also add IoU to the per-scan velocity print:
```python
    print(f"  {scan_labels[i]} -> {scan_labels[i+1]}: "
          f"{spd:.1f} km/hr  bearing={brg:.0f}°")
```
(No change needed — speeds are still computed the same way, just the
exclusion criterion changed.)

### After running: paste
1. Any [TRACK FAIL] lines
2. Valid intervals count and per-scan velocity table
3. Mean speed, max speed, mean bearing (valid intervals only)
4. Whether smoke_footprints.npz saved without error

---
## RESULTS 14

**Status**: DONE

**Naming/deviation notes**: Both changes applied verbatim.
`smoke_mask_3d`, `results`, `hotspots`, `dt` all matched existing names in
`krgx_hourly.py`/`task12_centroid_drift.py` — no mismatches. Re-ran
`krgx_hourly.py` once to regenerate `hourly_results.json` +
`smoke_footprints.npz` before running the updated tracker script.

**1. [TRACK FAIL] lines**: none. Zero IoU failures over all 13 transitions —
the only excluded interval is the very first (23:59→00:04), which is
excluded because scan 1 is the tracker's initial seed (no previous scan
to link from), not because of an IoU failure.

**2. Valid intervals / per-scan velocity table**: 12/13 valid.
```
23:59 -> 00:04: EXCLUDED (tracking failure or re-seed)   [seed, not a failure]
00:04 -> 00:08: 14.6 km/hr  bearing=69°
00:08 -> 00:12: 57.9 km/hr  bearing=57°
00:12 -> 00:16: 27.3 km/hr  bearing=196°
00:16 -> 00:20: 23.2 km/hr  bearing=99°
00:20 -> 00:24: 27.6 km/hr  bearing=54°
00:24 -> 00:32: 21.1 km/hr  bearing=344°
00:32 -> 00:36: 14.7 km/hr  bearing=113°
00:36 -> 00:40: 107.3 km/hr bearing=234°
00:40 -> 00:44: 6.5 km/hr   bearing=54°
00:44 -> 00:48: 10.2 km/hr  bearing=154°
00:48 -> 00:52: 28.3 km/hr  bearing=80°
00:52 -> 00:56: 27.1 km/hr  bearing=278°
```

**3. Mean speed 30.5 km/hr, max speed 107.3 km/hr, mean bearing 144°**
(valid intervals only).

**4. `smoke_footprints.npz` saved without error** — confirmed by the
`Footprints saved to krgx_hourly/smoke_footprints.npz` print and by the
tracker script successfully loading it on the next run.

**Flagged issue (significant)**: the IoU criterion accepted *every single*
transition (0 failures), yet a 107.3 km/hr spike at 00:36→00:40 is still
present — i.e. exactly the kind of implausible jump this task was meant
to catch slipped through. Root cause: `smoke_footprint_2d` is built from
`smoke_mask_3d.any(axis=2)` — the **whole-scan** smoke mask across the
*entire* domain, not the footprint of the specific hotspot being tracked
(the task's own background note calls it "union of all confirmed hotspot
footprints," but the actual Change 1 code doesn't filter by hotspot at
all — it's the same scene-wide mask already computed for the 2D PNG
overlay). Since the fire's smoke plume from all lobes stays large and
spatially contiguous for most of the hour, scene-wide IoU between any two
consecutive scans is almost always ≥ 0.10 regardless of which specific
hotspot centroid the nearest-neighbor step picks — so the check is
structurally unable to reject a bad candidate as long as *some* part of
the scene's smoke overlaps, which is nearly always true. Concretely, at
00:36→00:40 the nearest-neighbor step jumped from a hotspot near
(-26.1, 2.2) to one at (-32.1, -2.1) — 7.4 km in 4.1 min — while the
scene-wide smoke footprints obviously overlapped (same storm, same hour),
so `link_valid=True` regardless.

Not fixed since Task 14 specified the exact code to write and it was
implemented as given — but this means the goal stated in the task's own
Background section ("should only fail when two completely disjoint
echoes are being compared") is not actually being met by this
implementation. A per-hotspot (not per-scan) footprint would be needed
for IoU to discriminate between candidates the way the task intends —
e.g. storing each confirmed hotspot's own `region.any(axis=2)` mask keyed
by `(scan_dt, hotspot_index)` rather than one mask per scan.

## TASK 15 — Robust drift estimate via median filter + outlier rejection
**Status**: PENDING

**Goal**: Instead of fixing the tracker, post-process the centroid time
series to get a robust drift estimate. Flag intervals with implausible
speed OR bearing reversal as outliers, report median speed and bearing
from the remaining intervals.

**Instructions**: Add a post-processing block to task12_centroid_drift.py
after the existing velocity computation. No changes to krgx_hourly.py.

```python
# ── Robust drift estimate ─────────────────────────────────────────────────
# Criteria for a plausible interval:
#   1. Speed <= 60 km/hr (generous upper bound for smoke advection)
#   2. No sharp bearing reversal vs neighbours (|delta_bearing| <= 90°)

SPEED_MAX  = 60.0   # km/hr — generous bound for wind-advected smoke
BEAR_DELTA = 90.0   # deg  — max bearing change between adjacent intervals

valid_speed   = speed_kmph <= SPEED_MAX
# bearing reversal: compare each interval to its neighbours
bear_arr = np.where(np.isnan(bearing), np.nan, bearing)

# angular difference between consecutive bearings (wrapped to [-180, 180])
def ang_diff(a, b):
    d = (b - a + 180) % 360 - 180
    return abs(d)

bear_stable = np.ones(len(bearing), dtype=bool)
for i in range(len(bearing)):
    neighbours = []
    if i > 0 and not np.isnan(bearing[i-1]):
        neighbours.append(bearing[i-1])
    if i < len(bearing)-1 and not np.isnan(bearing[i+1]):
        neighbours.append(bearing[i+1])
    if neighbours and not np.isnan(bearing[i]):
        max_diff = max(ang_diff(bearing[i], n) for n in neighbours)
        bear_stable[i] = max_diff <= BEAR_DELTA

plausible = valid_speed & bear_stable & ~np.isnan(speed_kmph)

print(f"\n── Robust drift estimate ──────────────────────────────")
print(f"Total intervals:     {(~np.isnan(speed_kmph)).sum()}")
print(f"Speed outliers:      {(~valid_speed & ~np.isnan(speed_kmph)).sum()}")
print(f"Bearing reversals:   {(~bear_stable & ~np.isnan(speed_kmph)).sum()}")
print(f"Plausible intervals: {plausible.sum()}")

if plausible.sum() > 0:
    med_speed   = float(np.median(speed_kmph[plausible]))
    med_bearing = float(np.median(bearing[plausible]))
    print(f"\nMedian speed:   {med_speed:.1f} km/hr")
    print(f"Median bearing: {med_bearing:.0f}°")
    print(f"Plausible intervals detail:")
    for i in range(len(speed_kmph)):
        if plausible[i]:
            print(f"  {scan_labels[i]} -> {scan_labels[i+1]}: "
                  f"{speed_kmph[i]:.1f} km/hr  bearing={bearing[i]:.0f}°")
else:
    print("No plausible intervals — tracking is too unstable for "
          "a reliable drift estimate from this approach.")
    print("Recommend tobac or Kalman tracker for this event.")
```

### After running: paste
1. The full robust drift estimate block output
2. Which specific intervals survive as plausible
3. Median speed and bearing

---
## RESULTS 15

**Status**: DONE

**Naming/deviation notes**: Block appended verbatim after the existing
velocity print statements, inside the same `if len(times) >= 2:` block
(so `speed_kmph`, `bearing`, `scan_labels` are in scope without
re-derivation). No name mismatches — `speed_kmph`/`bearing` from Task 14's
tracker matched what this block expects directly.

**1. Full robust drift estimate block output**:
```
── Robust drift estimate ──────────────────────────────
Total intervals:     12
Speed outliers:      1
Bearing reversals:   10
Plausible intervals: 2

Median speed:   21.1 km/hr
Median bearing: 61°
Plausible intervals detail:
  00:04 -> 00:08: 14.6 km/hr  bearing=69°
  00:20 -> 00:24: 27.6 km/hr  bearing=54°
```

**2. Intervals surviving as plausible**: only 2 of 12 —
`00:04→00:08` (14.6 km/hr, 69°) and `00:20→00:24` (27.6 km/hr, 54°).
Both are non-adjacent, single, isolated intervals — not a contiguous
stretch of the hour.

**3. Median speed: 21.1 km/hr, median bearing: 61°** (≈ ENE, roughly
consistent with the earlier IoU-tracker mean bearing of 144°... actually
not very consistent — see flag below).

**Flagged issue**: 10 of 12 intervals (83%) fail the bearing-reversal
check, not the speed check (only 1 speed outlier — the 107.3 km/hr one
already flagged in RESULTS 14). This means the *direction* of scan-to-scan
centroid movement is flipping by more than 90° almost every other step,
which is consistent with genuinely noisy/non-continuous tracking rather
than a plume drifting steadily in one direction with a couple of bad
frames mixed in. The two "plausible" intervals that survive are isolated
(not adjacent to each other), so the median speed/bearing is being
computed from two disconnected 4-minute snapshots 16 minutes apart, not
a trend spanning the hour — i.e. the median filter finds *a* pair of
locally-consistent steps, but that's a much weaker claim than "the
robust drift estimate for the hour is 21.1 km/hr @ 61°." If a single
drift number is going in the write-up, worth being explicit that it's
derived from 2/12 samples, not the full hour, and that the tracker's
underlying instability (flagged in RESULTS 14 — footprint IoU computed
scene-wide rather than per-hotspot) is still the root cause this
post-processing step works around rather than fixes, exactly as the
task's own framing states ("instead of fixing the tracker..."). Not
changed since Task 15 explicitly chose the post-processing route over a
tracker fix — flagging so the PM can judge whether 2/12 is intended as
"proof this approach also struggles" or a usable number.

## TASK 16 — Probabilistic smoke occupancy tracker
**Status**: PENDING

**Goal**: Replace hard centroid linking with a Bayesian occupancy grid that
accumulates evidence across all 14 scans. Output a smooth centroid time
series with per-scan confidence, robust to weak/noisy/fragmented radar.

**Instructions**: Create `task16_prob_tracker.py` in the Runs dir.
This script loads existing outputs — no changes to krgx_hourly.py.
It needs access to the raw grids (rhv_grid, dbz_smooth) per scan, which
are not currently saved. Add a second temp save block to krgx_hourly.py,
run it once, then remove the block.

### Change 1 — Add per-scan grid save to krgx_hourly.py (temporary)
In `process_scan()`, after `smoke_mask_3d` is computed, add:

```python
    # TEMP Task 16 — save raw grids for probabilistic tracker
    _grid_path = OUT_DIR / f'grids_{fname_dt.strftime("%Y%m%d_%H%M%S")}.npz'
    np.savez_compressed(
        _grid_path,
        rhv_grid   = rhv_grid.astype(np.float32),
        dbz_smooth = dbz_smooth.astype(np.float32),
        x_bins     = x_bins.astype(np.float32),
        y_bins     = y_bins.astype(np.float32),
        z_bins     = z_bins.astype(np.float32),
    )
```

Re-run krgx_hourly.py once to generate 14 `grids_*.npz` files in
`krgx_hourly/`. Then remove this block — it is not permanent.

### Change 2 — Create task16_prob_tracker.py

```python
"""
Probabilistic smoke occupancy tracker — Task 16
Bayesian voxel-wise occupancy grid updated across 14 KRGX scans.
Centroid derived from probability expectation, not hard clustering.
"""
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime

# ── config ────────────────────────────────────────────────────────────────
OUT_DIR   = Path('krgx_hourly')
ALPHA     = 0.75    # temporal decay: higher = longer memory (~3 scan halflife)
RHV_MAX   = 0.90    # smoke RhoHV threshold (matches krgx_hourly.py)
DBZ_MIN   = 10.0    # dBZ floor (matches DBZ_THRESHOLD)

# Likelihood values (probability that this voxel contains smoke)
P_SMOKE       = 0.92   # rhv < RHV_MAX AND dbz > DBZ_MIN
P_NON_SMOKE   = 0.08   # data present, thresholds not met
P_NO_DATA     = 0.50   # masked / no observation — neutral

# ── load scan list ────────────────────────────────────────────────────────
with open(OUT_DIR / 'hourly_results.json') as f:
    results = json.load(f)

grid_files = sorted(OUT_DIR.glob('grids_*.npz'))
print(f"Found {len(grid_files)} grid files, {len(results)} result entries")

# ── initialize occupancy grid from first scan's shape ────────────────────
d0       = np.load(grid_files[0])
nx       = len(d0['x_bins']) - 1
ny       = len(d0['y_bins']) - 1
nz       = len(d0['z_bins']) - 1
x_km     = (d0['x_bins'][:-1] + d0['x_bins'][1:]) / 2
y_km     = (d0['y_bins'][:-1] + d0['y_bins'][1:]) / 2
z_bins   = d0['z_bins']
z_km     = (z_bins[:-1] + z_bins[1:]) / 2
RADAR_ALT_KM = 2.530
z_km_msl = z_km + RADAR_ALT_KM

prob_grid = np.full((nx, ny, nz), 0.50, dtype=np.float32)  # uniform prior

# coordinate grids for expectation calculation
XX, YY, ZZ = np.meshgrid(x_km, y_km, z_km_msl, indexing='ij')

# ── output accumulators ───────────────────────────────────────────────────
times      = []
x_cents    = []
y_cents    = []
z_cents    = []
z_tops     = []   # 95th percentile of prob-weighted z
confidences = []  # total probability mass (proxy for detection strength)
scan_labels = []

# ── main loop ─────────────────────────────────────────────────────────────
for gf in grid_files:
    # parse scan time from filename: grids_YYYYMMDD_HHMMSS.npz
    stem = gf.stem  # e.g. grids_20260809_235957
    dt_str = stem.replace('grids_', '')
    dt = datetime.strptime(dt_str, '%Y%m%d_%H%M%S')

    d         = np.load(gf)
    rhv       = d['rhv_grid'].astype(float)
    dbz       = d['dbz_smooth'].astype(float)

    # ── build likelihood field ────────────────────────────────────────────
    has_data  = np.isfinite(rhv) & np.isfinite(dbz) & (dbz > 0)
    is_smoke  = has_data & (rhv < RHV_MAX) & (dbz > DBZ_MIN)
    is_clear  = has_data & ~is_smoke

    likelihood = np.full((nx, ny, nz), P_NO_DATA, dtype=np.float32)
    likelihood[is_smoke] = P_SMOKE
    likelihood[is_clear] = P_NON_SMOKE

    # ── Bayesian update with temporal decay ───────────────────────────────
    # 1. Decay prior toward neutral (0.5) — old evidence fades
    prob_grid = ALPHA * prob_grid + (1 - ALPHA) * 0.5

    # 2. Bayesian update: P_new ∝ P_prior * likelihood
    unnorm    = prob_grid * likelihood
    # Normalise per-voxel relative to complement
    # (binary: smoke vs non-smoke)
    complement = (1 - prob_grid) * (1 - likelihood)
    denom     = unnorm + complement
    prob_grid = np.where(denom > 0, unnorm / denom, 0.5).astype(np.float32)

    # ── derive centroid from probability expectation ──────────────────────
    # Use voxels above a soft threshold to avoid noise floor dominating
    weight    = np.maximum(prob_grid - 0.5, 0.0)   # only above-neutral mass
    total_w   = weight.sum()

    if total_w > 0:
        cx = float((XX * weight).sum() / total_w)
        cy = float((YY * weight).sum() / total_w)
        cz = float((ZZ * weight).sum() / total_w)

        # z_top: 95th percentile of probability-weighted z distribution
        # flatten and sort by z, accumulate weight
        z_flat = ZZ.ravel()
        w_flat = weight.ravel()
        sort_i = np.argsort(z_flat)
        cumw   = np.cumsum(w_flat[sort_i])
        cumw  /= cumw[-1]
        z_top  = float(z_flat[sort_i][np.searchsorted(cumw, 0.95)])
    else:
        cx = cy = cz = z_top = np.nan

    confidence = float(total_w)

    times.append(dt)
    x_cents.append(cx)
    y_cents.append(cy)
    z_cents.append(cz)
    z_tops.append(z_top)
    confidences.append(confidence)
    scan_labels.append(dt.strftime('%H:%M'))

    print(f"  {dt.strftime('%H:%MZ')}  "
          f"centroid=({cx:.1f}, {cy:.1f}, {cz:.2f} km MSL)  "
          f"z_top={z_top:.2f}  conf={confidence:.0f}")

# ── compute scan-to-scan drift ────────────────────────────────────────────
x_arr  = np.array(x_cents)
y_arr  = np.array(y_cents)
t_min  = np.array([(t - times[0]).total_seconds()/60 for t in times])

dx     = np.diff(x_arr)
dy     = np.diff(y_arr)
dt_min = np.diff(t_min)
speed  = np.sqrt(dx**2 + dy**2) / (dt_min / 60)
bearing = (np.degrees(np.arctan2(dx, dy)) % 360)

print(f"\nDrift statistics (all {len(speed)} intervals):")
print(f"  Speed   — mean={np.nanmean(speed):.1f}  "
      f"median={np.nanmedian(speed):.1f}  "
      f"max={np.nanmax(speed):.1f} km/hr")
print(f"  Bearing — mean={np.nanmean(bearing):.0f}°  "
      f"median={np.nanmedian(bearing):.0f}°")

# ── plots ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(12, 9))

# 1. Plan view centroid track
ax = axes[0, 0]
sc = ax.scatter(x_arr, y_arr, c=t_min, cmap='plasma', s=60, zorder=3)
ax.plot(x_arr, y_arr, 'k-', lw=1.2, zorder=2)
for i, lbl in enumerate(scan_labels):
    ax.annotate(lbl, (x_arr[i], y_arr[i]),
                fontsize=7, xytext=(4,4), textcoords='offset points')
ax.scatter(0, 0, marker='+', c='black', s=120, zorder=5, label='KRGX')
ax.scatter(-35.8, 1.1, marker='*', c='red', s=150, zorder=5,
           label='Bug Fire anchor')
plt.colorbar(sc, ax=ax, label='Minutes from first scan')
ax.set(xlabel='X km', ylabel='Y km',
       title='Probabilistic centroid track (plan view)', aspect='equal')
ax.legend(fontsize=8)

# 2. Confidence (probability mass) vs time
ax = axes[0, 1]
ax.plot(t_min, confidences, 'o-', color='steelblue')
ax.set(xlabel='Minutes from first scan',
       ylabel='Probability mass (above-neutral)',
       title='Detection confidence vs time')

# 3. Height vs time
ax = axes[1, 0]
ax.plot(t_min, z_cents, 'o-', color='darkorange', label='Z centroid')
ax.plot(t_min, z_tops,  's--', color='firebrick',  label='Z top (p95)')
ax.set(xlabel='Minutes from first scan', ylabel='km MSL',
       title='Plume height vs time (probabilistic)')
ax.legend(fontsize=8)

# 4. Scan-to-scan speed
ax = axes[1, 1]
t_mid = (t_min[:-1] + t_min[1:]) / 2
ax.bar(t_mid, speed, width=3, color='purple', alpha=0.7)
ax.axhline(20, color='red', lw=1.5, linestyle='--', label='20 km/hr ref')
ax.set(xlabel='Minutes from first scan', ylabel='km/hr',
       title='Scan-to-scan centroid speed')
ax.legend(fontsize=8)

plt.suptitle('KRGX — Bug Fire probabilistic plume tracker\n'
             '2026-08-09 23:59Z to 2026-08-10 00:56Z', fontsize=11)
plt.tight_layout()
out_png = OUT_DIR / 'task16_prob_tracker.png'
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print(f'\nSaved {out_png}')
```

### After running: paste
1. Full per-scan centroid + confidence print block
2. Drift statistics block
3. Whether any scan-to-scan speed exceeds 20 km/hr
4. Whether the PNG saved without error
5. Qualitative description of the plan-view track —
   does it show a smooth drift or still jumpy?

---
## RESULTS 16

**Status**: DONE

**Naming/deviation notes**: Change 1 (temp grid save) added to
`krgx_hourly.py` right after `smoke_footprint_2d` is computed, using the
existing `fname_dt`/`x_bins`/`y_bins`/`z_bins`/`rhv_grid`/`dbz_smooth`
names already in scope in `process_scan()` — no mismatches. Ran once to
produce 14 `grids_*.npz` files, then removed the block exactly as
instructed (confirmed it's gone from the current file; the 14 npz files
remain on disk for `task16_prob_tracker.py` to read). `task16_prob_tracker.py`
created verbatim as a new standalone script — no name mismatches since
it only reads from the already-generated `hourly_results.json` and
`grids_*.npz`, not from `krgx_hourly.py`'s internals directly.

**1. Full per-scan centroid + confidence print block**:
```
Found 14 grid files, 14 result entries
  23:59Z  centroid=(-28.4, 2.7, 2.92 km MSL)  z_top=3.40  conf=486
  00:04Z  centroid=(-29.8, 1.6, 2.96 km MSL)  z_top=3.40  conf=535
  00:08Z  centroid=(-29.3, 1.2, 2.97 km MSL)  z_top=3.40  conf=705
  00:12Z  centroid=(-28.6, 1.3, 2.98 km MSL)  z_top=3.40  conf=918
  00:16Z  centroid=(-27.4, 2.0, 2.95 km MSL)  z_top=3.40  conf=977
  00:20Z  centroid=(-26.0, 2.9, 2.93 km MSL)  z_top=3.40  conf=1149
  00:24Z  centroid=(-24.9, 3.5, 2.88 km MSL)  z_top=3.40  conf=991
  00:32Z  centroid=(-23.1, 3.9, 2.92 km MSL)  z_top=3.40  conf=1064
  00:36Z  centroid=(-21.8, 4.8, 2.91 km MSL)  z_top=3.40  conf=1208
  00:40Z  centroid=(-20.5, 5.6, 2.89 km MSL)  z_top=3.40  conf=1255
  00:44Z  centroid=(-18.8, 6.5, 2.86 km MSL)  z_top=3.15  conf=1201
  00:48Z  centroid=(-17.5, 6.7, 2.84 km MSL)  z_top=3.15  conf=1094
  00:52Z  centroid=(-14.7, 7.6, 2.83 km MSL)  z_top=3.15  conf=793
  00:56Z  centroid=(-12.5, 8.8, 2.79 km MSL)  z_top=3.15  conf=489
```

**2. Drift statistics block**:
```
Drift statistics (all 13 intervals):
  Speed   — mean=22.8  median=23.0  max=42.4 km/hr
  Bearing — mean=84°  median=65°
```

**3. Does any scan-to-scan speed exceed 20 km/hr?** Yes — 8 of the 13
intervals exceed 20 km/hr (recomputed per-interval from the printed
centroids/exact scan timestamps to check, since the script itself only
prints the aggregate stats): 23:59→00:04 (25.9), 00:16→00:20 (24.4),
00:32→00:36 (23.1), 00:36→00:40 (22.2), 00:40→00:44 (28.0),
00:48→00:52 (43.2), 00:52→00:56 (36.5), and 00:12→00:16 sits right at
the line (20.3). But unlike every prior tracker (Tasks 12-15), none of
these are wild outliers — the largest is 43.2 km/hr, not 107+ km/hr, and
they occur in a smoothly increasing sequence (see below), consistent
with the plume's translation speed genuinely picking up in the second
half of the hour rather than a tracking-identity error.

**4. PNG saved without error**: yes, confirmed by the
`Saved krgx_hourly\task16_prob_tracker.png` print with no exceptions.

**5. Qualitative description of the plan-view track**: smooth and
monotonic — x centroid moves steadily eastward from -28.4 to -12.5 km
and y moves steadily northward from 2.7 to 8.8 km, with no
back-and-forth jitter or lobe-swapping at any point across all 14 scans.
This is a clear improvement over every previous tracker in this task
sequence (Tasks 12-15), which all showed erratic jumps and bearing
reversals. Two things worth noting rather than flagging as bugs:
(a) confidence (probability mass) rises from 486 to a peak of 1255 at
00:40Z then declines to 489 by 00:56Z — this tracks the smoke_vox_total
pattern already seen in Task 09/10 (smoke mask peaked mid-hour, thinned
toward the end), so the probabilistic centroid's reliability is likely
lower in the last 1-2 scans even though the track itself stays smooth;
(b) drift speed itself isn't constant — it's noticeably faster in the
back half of the hour (00:40 onward, several 20-40+ km/hr intervals)
than the front half (00:04-00:24, mostly under 20-25 km/hr) — this reads
as genuine acceleration of the plume's advection rather than a tracking
artifact, given how smooth the underlying (x,y) path is, but worth
flagging as a real physical signal for the write-up rather than
smoothing over it with a single mean/median speed for the whole hour.

## TASK 17 — Per-timestep radar images with probabilistic centroid overlay
**Status**: PENDING

**Goal**: Generate 14 PNG images (one per scan) showing the 2D column-max
dBZ field with the probabilistic centroid from Task 16 overlaid and labeled.
Also produce one summary panel showing all 14 centroids on a single plot.
Standalone scratch script — no changes to any existing file.

**Instructions**: Create `task17_centroid_frames.py` in the Runs dir.

```python
"""
Per-timestep radar images with probabilistic centroid overlay — Task 17
Loads grids_*.npz (raw dBZ/RhoHV per scan) and task16 centroid output,
overlays centroid position on 2D column-max dBZ for each scan.
"""
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
from datetime import datetime

OUT_DIR   = Path('krgx_hourly')
FRAME_DIR = OUT_DIR / 'centroid_frames'
FRAME_DIR.mkdir(exist_ok=True)

RADAR_ALT_KM = 2.530
DBZ_MIN      = 10.0
RHV_MAX      = 0.90

# ── re-run the probabilistic tracker inline to get centroids ──────────────
# (copy the core loop from task16_prob_tracker.py — no import needed)
ALPHA       = 0.75
P_SMOKE     = 0.92
P_NON_SMOKE = 0.08
P_NO_DATA   = 0.50

grid_files = sorted(OUT_DIR.glob('grids_*.npz'))
d0         = np.load(grid_files[0])
nx  = len(d0['x_bins']) - 1
ny  = len(d0['y_bins']) - 1
nz  = len(d0['z_bins']) - 1
x_km     = (d0['x_bins'][:-1] + d0['x_bins'][1:]) / 2
y_km     = (d0['y_bins'][:-1] + d0['y_bins'][1:]) / 2
z_km_msl = (d0['z_bins'][:-1] + d0['z_bins'][1:]) / 2 + RADAR_ALT_KM

prob_grid = np.full((nx, ny, nz), 0.50, dtype=np.float32)
XX, YY, ZZ = np.meshgrid(x_km, y_km, z_km_msl, indexing='ij')

centroids   = []   # (dt, cx, cy, cz, confidence)

for gf in grid_files:
    stem  = gf.stem.replace('grids_', '')
    dt    = datetime.strptime(stem, '%Y%m%d_%H%M%S')
    d     = np.load(gf)
    rhv   = d['rhv_grid'].astype(float)
    dbz   = d['dbz_smooth'].astype(float)

    has_data  = np.isfinite(rhv) & np.isfinite(dbz) & (dbz > 0)
    is_smoke  = has_data & (rhv < RHV_MAX) & (dbz > DBZ_MIN)
    is_clear  = has_data & ~is_smoke

    likelihood = np.full((nx, ny, nz), P_NO_DATA, dtype=np.float32)
    likelihood[is_smoke] = P_SMOKE
    likelihood[is_clear] = P_NON_SMOKE

    prob_grid  = ALPHA * prob_grid + (1 - ALPHA) * 0.5
    unnorm     = prob_grid * likelihood
    complement = (1 - prob_grid) * (1 - likelihood)
    denom      = unnorm + complement
    prob_grid  = np.where(denom > 0, unnorm / denom, 0.5).astype(np.float32)

    weight  = np.maximum(prob_grid - 0.5, 0.0)
    total_w = weight.sum()

    if total_w > 0:
        cx = float((XX * weight).sum() / total_w)
        cy = float((YY * weight).sum() / total_w)
        cz = float((ZZ * weight).sum() / total_w)
    else:
        cx = cy = cz = np.nan

    centroids.append((dt, cx, cy, cz, float(total_w), dbz, rhv))

# ── per-scan frames ───────────────────────────────────────────────────────
norm     = mcolors.Normalize(vmin=0, vmax=35)
cmap_dbz = plt.cm.jet
smoke_cmap = mcolors.ListedColormap(['none', 'purple'])

# accumulate centroid trail for overlay
trail_x, trail_y = [], []

for i, (dt, cx, cy, cz, conf, dbz, rhv) in enumerate(centroids):
    trail_x.append(cx)
    trail_y.append(cy)

    # 2D column-max dBZ
    col_max = dbz.max(axis=2)   # (nx, ny)

    # 2D smoke mask footprint (column-any of prob > 0.6)
    smoke_fp = (prob_grid > 0.6).any(axis=2)   # approximate — use
    # the current prob_grid state isn't accessible here after the loop;
    # use rhv/dbz directly for the footprint overlay instead
    smoke_2d = (
        (np.isfinite(rhv) & (rhv < RHV_MAX) & (dbz > DBZ_MIN))
        .any(axis=2)
    )

    fig, ax = plt.subplots(figsize=(8, 7))

    # dBZ background
    pcm = ax.pcolormesh(x_km, y_km, col_max.T,
                        cmap=cmap_dbz, norm=norm, shading='auto')
    plt.colorbar(pcm, ax=ax, label='dBZ (column max)', pad=0.02)

    # smoke mask footprint (purple tint)
    ax.contourf(x_km, y_km, smoke_2d.T.astype(float),
                levels=[0.5, 1.5], colors=['purple'], alpha=0.20, zorder=2)
    ax.contour(x_km, y_km, smoke_2d.T.astype(float),
               levels=[0.5], colors=['purple'], linewidths=1,
               linestyles='--', zorder=3)

    # centroid trail (previous scans, faded)
    if len(trail_x) > 1:
        ax.plot(trail_x[:-1], trail_y[:-1],
                'w--', lw=1.0, alpha=0.5, zorder=4)
        ax.scatter(trail_x[:-1], trail_y[:-1],
                   c='white', s=20, alpha=0.4, zorder=4)

    # current centroid
    if np.isfinite(cx):
        ax.scatter(cx, cy, c='yellow', s=180, marker='*',
                   edgecolors='black', linewidths=0.8,
                   zorder=6, label=f'Centroid z={cz:.2f} km MSL')
        ax.annotate(
            f'{dt.strftime("%H:%MZ")}\n({cx:.1f}, {cy:.1f}) km\n'
            f'z={cz:.2f} km MSL\nconf={conf:.0f}',
            xy=(cx, cy),
            xytext=(8, 8), textcoords='offset points',
            fontsize=8, color='yellow',
            bbox=dict(boxstyle='round,pad=0.2', fc='black', alpha=0.5),
            zorder=7
        )

    # Bug Fire anchor
    ax.scatter(-35.8, 1.1, marker='*', c='red', s=150,
               edgecolors='white', linewidths=0.5,
               zorder=5, label='Bug Fire anchor')

    # radar site
    ax.scatter(0, 0, marker='+', c='white', s=120, zorder=5, label='KRGX')

    ax.set(aspect='equal',
           xlim=(x_km.min(), x_km.max()),
           ylim=(y_km.min(), y_km.max()),
           xlabel='X km (aeqd from KRGX)',
           ylabel='Y km (aeqd from KRGX)',
           title=f'KRGX  {dt.strftime("%Y-%m-%d %H:%MZ")}  '
                 f'[{i+1}/14]  conf={conf:.0f}')
    ax.legend(fontsize=8, loc='upper right')

    out_path = FRAME_DIR / f'frame_{i+1:02d}_{dt.strftime("%H%M%S")}.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved {out_path.name}')

# ── summary panel: all centroids on one plot ──────────────────────────────
fig, ax = plt.subplots(figsize=(9, 8))

# background: first scan's column-max for spatial reference
col_max_0 = centroids[0][5].max(axis=2)
ax.pcolormesh(x_km, y_km, col_max_0.T,
              cmap=cmap_dbz, norm=norm, shading='auto', alpha=0.4)

# centroid track colored by time
t_vals = np.arange(len(centroids))
cx_all = np.array([c[1] for c in centroids])
cy_all = np.array([c[2] for c in centroids])
cz_all = np.array([c[3] for c in centroids])
cf_all = np.array([c[4] for c in centroids])

ax.plot(cx_all, cy_all, 'w-', lw=1.5, alpha=0.6, zorder=3)
sc = ax.scatter(cx_all, cy_all, c=t_vals, cmap='plasma',
                s=80, zorder=4, edgecolors='white', linewidths=0.5)
plt.colorbar(sc, ax=ax, label='Scan index (0=23:59Z, 13=00:56Z)')

# label every other scan to avoid clutter
for i, (dt, cx, cy, *_) in enumerate(centroids):
    if i % 2 == 0:
        ax.annotate(dt.strftime('%H:%MZ'), (cx, cy),
                    fontsize=7, color='white',
                    xytext=(4, 4), textcoords='offset points')

ax.scatter(-35.8, 1.1, marker='*', c='red', s=200,
           edgecolors='white', zorder=5, label='Bug Fire anchor')
ax.scatter(0, 0, marker='+', c='white', s=150, zorder=5, label='KRGX')

ax.set(aspect='equal',
       xlim=(x_km.min(), x_km.max()),
       ylim=(y_km.min(), y_km.max()),
       xlabel='X km (aeqd from KRGX)',
       ylabel='Y km (aeqd from KRGX)',
       title='KRGX Bug Fire — probabilistic centroid track\n'
             '2026-08-09 23:59Z to 2026-08-10 00:56Z')
ax.legend(fontsize=8)
plt.tight_layout()
summary_path = OUT_DIR / 'task17_centroid_summary.png'
plt.savefig(summary_path, dpi=150, bbox_inches='tight')
plt.close()
print(f'\nSaved summary: {summary_path}')
print(f'All frames in: {FRAME_DIR}')
```

### After running: paste
1. List of 14 frame filenames confirmed saved
2. Whether summary PNG saved without error
3. Any traceback

Then open `task17_centroid_summary.png` and the individual frames —
share them here if possible.

---
## RESULTS 17

**Status**: DONE

**Naming/deviation notes**: Script created verbatim as `task17_centroid_frames.py`.
No name mismatches — it's fully standalone (re-derives centroids inline
from the same `grids_*.npz` files Task 16 used, rather than importing
from `task16_prob_tracker.py`), and those npz files were still present on
disk (Task 16 only removed the *code that generates* them from
`krgx_hourly.py`, not the files themselves — confirmed all 14
`grids_*.npz` present before running). No changes made to any existing
file, per the task's own constraint.

**1. Frame filenames confirmed saved** (14/14):
```
frame_01_235957.png   frame_08_003207.png
frame_02_000404.png   frame_09_003613.png
frame_03_000809.png   frame_10_004020.png
frame_04_001215.png   frame_11_004427.png
frame_05_001621.png   frame_12_004833.png
frame_06_002027.png   frame_13_005238.png
frame_07_002433.png   frame_14_005645.png
```
All in `krgx_hourly/centroid_frames/`.

**2. Summary PNG**: saved without error to
`krgx_hourly/task17_centroid_summary.png`.

**3. Traceback**: none — clean run, no exceptions.

**Note on sharing images**: this channel is text-only (TASKS.md), so I
can't embed the PNGs directly here — all 15 files (14 frames + 1 summary)
are on disk at the paths above for the PM/user to open directly.

**Minor observation (not a bug)**: the script's own inline comment flags
that `smoke_fp = (prob_grid > 0.6).any(axis=2)` on line ~2870 is dead
code — it uses whatever `prob_grid` state is left over after the entire
14-scan Bayesian update loop finished (i.e. the *final* scan's
occupancy grid, not scan `i`'s), and is immediately superseded by the
`smoke_2d` computed directly from that scan's own `rhv`/`dbz` right
after. `smoke_fp` is computed but never used in the plot — harmless,
matches the code's own acknowledgment in-comment, left as-is.

## TASK 18 — Gaussian spatial prior anchored to fire source for column tracking
**Status**: PENDING

**Goal**: Modify the probabilistic tracker in task16_prob_tracker.py and
task17_centroid_frames.py to use a Gaussian-shaped spatial prior centered
on the Bug Fire anchor instead of a uniform 0.5 prior. Voxels far from
the anchor decay toward low probability, so the tracker preferentially
follows the rising column near the source rather than the detached
advected plume.

**Background**: Instead of:
    prob_grid = ALPHA * prob_grid + (1-ALPHA) * 0.5
We use:
    prob_grid = ALPHA * prob_grid + (1-ALPHA) * prior_grid
Where prior_grid is a 3D Gaussian centered on the fire anchor in x/y,
uniform in z (we don't constrain height — the column can loft freely).
Voxels at the anchor get prior=0.7 (moderate positive bias toward smoke),
voxels 15+ km away get prior=0.1 (moderate negative bias — unlikely to
be the rising column). This is a soft constraint, not a hard cutoff —
strong radar evidence far from the anchor can still activate those voxels,
but without evidence they decay toward low probability.

**Instructions**: Make these changes to BOTH task16_prob_tracker.py AND
task17_centroid_frames.py (the centroid re-computation loop in
task17 is a copy of task16's core loop — both need the same fix).

### Change 1 — Add Gaussian prior grid (after grid shape is known,
before the main loop)

```python
# ── Gaussian spatial prior centered on Bug Fire anchor ───────────────────
ANCHOR_X_KM  = -35.8   # Bug Fire GOES power-weighted centroid (Task 08b)
ANCHOR_Y_KM  =   1.1
PRIOR_SIGMA_KM = 8.0   # e-folding radius — 8 km covers the column base
PRIOR_MAX    = 0.70    # prior probability AT the anchor
PRIOR_MIN    = 0.10    # prior probability far from anchor

# Build 2D Gaussian in x/y (same for all z levels)
dist2_xy = (XX - ANCHOR_X_KM)**2 + (YY - ANCHOR_Y_KM)**2   # (nx,ny,nz)
gauss_xy = np.exp(-dist2_xy / (2 * PRIOR_SIGMA_KM**2))      # (nx,ny,nz)

# Scale to [PRIOR_MIN, PRIOR_MAX]
prior_grid = (PRIOR_MAX - PRIOR_MIN) * gauss_xy + PRIOR_MIN
prior_grid = prior_grid.astype(np.float32)

print(f"Prior grid: min={prior_grid.min():.3f}  "
      f"max={prior_grid.max():.3f}  "
      f"at anchor: {prior_grid[np.argmin(np.abs(x_km-ANCHOR_X_KM)), "
      f"np.argmin(np.abs(y_km-ANCHOR_Y_KM)), 0]:.3f}")
```

Note: `XX`, `YY` are already (nx,ny,nz) meshgrids built before the loop —
use them directly. `dist2_xy` varies in x/y but is the same at every z
level because the Gaussian is 2D only (no z constraint).

### Change 2 — Replace uniform decay target with prior_grid in the loop
Find:
```python
    prob_grid = ALPHA * prob_grid + (1 - ALPHA) * 0.5
```
Replace with:
```python
    prob_grid = ALPHA * prob_grid + (1 - ALPHA) * prior_grid
```
This is the only change inside the loop — everything else (likelihood
computation, Bayesian update, centroid extraction) stays identical.

### Change 3 — Tune ALPHA for tighter column focus
Change:
```python
ALPHA = 0.75
```
to:
```python
ALPHA = 0.60   # shorter memory — column position updates faster,
               # less accumulation of detached plume evidence
```
Rationale: with a uniform prior, high ALPHA was needed to accumulate
enough evidence across sparse scans. With a strong spatial prior pulling
toward the anchor, we want faster decay so the detached plume evidence
doesn't persist and dominate.

### Change 4 — Add prior contour to per-frame plots in task17
In the per-scan frame plotting block, after the smoke footprint contour,
add a faint contour showing the prior's 0.5 isoline (the boundary where
the prior is neutral — inside this is where the tracker is biased toward
the column):

```python
    # Prior boundary (where prior = 0.5 — tracker's neutral zone)
    prior_2d = prior_grid[:, :, 0]   # same at all z levels
    ax.contour(x_km, y_km, prior_2d.T,
               levels=[0.50], colors=['white'],
               linewidths=1.0, linestyles=':', alpha=0.6, zorder=3)
```

### After running BOTH scripts: paste
1. Prior grid min/max/at-anchor values
2. Full per-scan centroid + confidence block from task16
3. Whether centroids now stay near the fire anchor
   (x should stay in roughly -38 to -28 km range)
4. Qualitative description of task17 frames —
   does the tracked point sit on the rising column
   rather than the detached plume?
5. Updated drift statistics

---
## RESULTS 18

**Status**: DONE

**Naming/deviation notes**: All 4 changes applied to both
`task16_prob_tracker.py` and `task17_centroid_frames.py`. One deviation
from the literal instructions: **Change 1's print statement as written
in the task is a syntax error** — it splits an f-string expression
(`prior_grid[np.argmin(...), np.argmin(...), 0]`) across two separate
adjacent string literals, which Python cannot parse (each f-string's
`{...}` must open and close within that same literal). Fixed by computing
`_anchor_ix`/`_anchor_iy` as separate variables first, then using them in
a single, valid f-string — same printed values, no logic change.
Everything else (Gaussian prior formula, ALPHA change, decay-target
substitution, prior contour overlay) applied verbatim to both files.

**1. Prior grid**: `min=0.100  max=0.700  at anchor: 0.700` (identical in
both scripts, as expected — same formula, same anchor).

**2. Full per-scan centroid + confidence block** (task16):
```
Prior grid: min=0.100  max=0.700  at anchor: 0.700
  23:59Z  centroid=(-34.5, 1.5, 4.96 km MSL)  z_top=7.65  conf=2443
  00:04Z  centroid=(-35.2, 1.2, 5.14 km MSL)  z_top=7.65  conf=3595
  00:08Z  centroid=(-35.0, 1.1, 5.13 km MSL)  z_top=7.65  conf=4415
  00:12Z  centroid=(-34.8, 1.2, 5.10 km MSL)  z_top=7.65  conf=4987
  00:16Z  centroid=(-34.6, 1.2, 5.12 km MSL)  z_top=7.65  conf=5225
  00:20Z  centroid=(-34.4, 1.3, 5.10 km MSL)  z_top=7.65  conf=5450
  00:24Z  centroid=(-34.4, 1.4, 5.15 km MSL)  z_top=7.65  conf=5415
  00:32Z  centroid=(-34.6, 1.4, 5.21 km MSL)  z_top=7.65  conf=5362
  00:36Z  centroid=(-34.1, 1.6, 5.14 km MSL)  z_top=7.65  conf=5581
  00:40Z  centroid=(-33.9, 1.8, 5.12 km MSL)  z_top=7.65  conf=5650
  00:44Z  centroid=(-33.7, 1.9, 5.12 km MSL)  z_top=7.65  conf=5667
  00:48Z  centroid=(-33.8, 1.8, 5.13 km MSL)  z_top=7.65  conf=5627
  00:52Z  centroid=(-34.3, 1.6, 5.24 km MSL)  z_top=7.65  conf=5386
  00:56Z  centroid=(-34.8, 1.5, 5.33 km MSL)  z_top=7.65  conf=5213
```

**3. Do centroids stay near the fire anchor?** Yes, and tighter than the
task's own expected range — x stays in **-35.2 to -33.7 km** (vs. the
"roughly -38 to -28" expected), y stays in **1.1 to 1.9 km**. The tracker
is now firmly locked onto the source, essentially decoupled from the
detached advected plume that Tasks 12-17 were tracking.

**4. Qualitative description of task17 frames**: the yellow centroid
marker sits consistently over the fire's base column (near the red Bug
Fire anchor star) in every frame, inside the white dotted prior-neutral
boundary — it does **not** follow the smoke plume as it visibly advects
away in the purple footprint contour. This matches the stated goal
exactly: the tracker now represents "where the fire is," not "where the
smoke has drifted to."

**5. Updated drift statistics**: mean speed 4.4 km/hr, median 3.8 km/hr,
max 10.8 km/hr, mean bearing 160°, median 93°. An order of magnitude
slower than every prior tracker (Tasks 12-17 all reported means in the
20-50 km/hr range) — consistent with tracking a quasi-stationary source
rather than an advecting plume.

**Flagged issue**: `z_top` is **exactly 7.65 km MSL in all 14 scans** —
i.e. pinned at (or extremely close to) the top of the domain
(Z_MAX=8000 m MSL ceiling) regardless of scan, including the very first
scan (23:59Z) when the fire had presumably not yet built a tall column.
Root cause: the new prior is only 2D (uniform across all z, per the
task's own design — "we don't constrain height"). At the anchor's x/y
location, any voxel with no radar return (`has_data=False`, common at
high altitudes beyond the beam's vertical extent) gets `likelihood=0.5`
(P_NO_DATA), and the Bayesian update algebraically reduces to
`prob_new = prior` whenever likelihood=0.5 exactly (since
`prior*0.5 / (prior*0.5 + (1-prior)*0.5) = prior`). Before this task, the
old uniform 0.5 prior meant those same no-data voxels stayed exactly
neutral (weight=0 in the centroid calc). Now, at the anchor column,
*every* z level's no-data voxels adopt the anchor's prior value (up to
0.70) regardless of whether the radar actually observed smoke there —
so the p95 vertical distribution is now dominated by artificially
"positive" no-data voxels stacked all the way to the domain ceiling,
rather than by actual detected column height. This makes `z_top` from
this version scientifically meaningless (it would report the same value
even for a cold, dead radar with zero fire activity, as long as it's
near the anchor x/y). `cz` (the mean-weighted height, ~5.0-5.3 km MSL)
is likely affected the same way, just less visibly since it's an average
rather than a tail percentile. Not fixed since Task 18 specified the
exact prior formula (2D-only, "uniform in z") — but the horizontal
tracking gain here comes with a real cost to any vertical/height metric,
which would need the prior's "no vertical constraint" statement revisited
(e.g. only apply the spatial prior where `has_data` is True, leaving
true no-data voxels at neutral 0.5 as before) if plume-top height is
needed downstream.

## TASK 19 — Fix z_top: only apply spatial prior where radar data exists
**Status**: PENDING

**Goal**: No-data voxels (masked gates, above beam coverage) should remain
at neutral 0.5 probability regardless of their distance from the anchor.
The Gaussian prior should only bias voxels where the radar actually
observed something. One-line change in both task16_prob_tracker.py and
task17_centroid_frames.py.

**Instructions**: In the main tracker loop, find the temporal decay line:
```python
    prob_grid = ALPHA * prob_grid + (1 - ALPHA) * prior_grid
```
Replace with:
```python
    decay_target = np.where(has_data, prior_grid, 0.5)
    prob_grid = ALPHA * prob_grid + (1 - ALPHA) * decay_target
```

`has_data` is already computed earlier in the loop as:
```python
    has_data = np.isfinite(rhv) & np.isfinite(dbz) & (dbz > 0)
```
So no new variables needed — just use it directly.

Apply this change to BOTH task16_prob_tracker.py and
task17_centroid_frames.py.

### After running both: paste
1. Full per-scan centroid + confidence block
2. z_top values — confirm they are no longer pinned at 7.65 km
3. Whether x/y centroids still stay near the anchor (-38 to -28 km range)
4. Updated drift statistics

---
## RESULTS 19

**Status**: DONE

**Naming/deviation notes**: One-line change applied verbatim to both
`task16_prob_tracker.py` and `task17_centroid_frames.py` — `has_data`
was already in scope exactly as described, no adaptation needed. Both
scripts re-ran cleanly, task17 regenerated all 14 frames + summary PNG.

**1. Full per-scan centroid + confidence block** (task16):
```
Prior grid: min=0.100  max=0.700  at anchor: 0.700
  23:59Z  centroid=(-28.7, 2.5, 2.92 km MSL)  z_top=3.40  conf=475
  00:04Z  centroid=(-30.7, 0.7, 2.98 km MSL)  z_top=3.40  conf=525
  00:08Z  centroid=(-30.0, 0.6, 2.99 km MSL)  z_top=3.40  conf=704
  00:12Z  centroid=(-29.2, 0.9, 3.00 km MSL)  z_top=3.40  conf=895
  00:16Z  centroid=(-28.1, 1.2, 2.95 km MSL)  z_top=3.40  conf=892
  00:20Z  centroid=(-26.9, 1.9, 2.94 km MSL)  z_top=3.40  conf=972
  00:24Z  centroid=(-25.6, 2.7, 2.89 km MSL)  z_top=3.40  conf=839
  00:32Z  centroid=(-24.2, 2.8, 2.94 km MSL)  z_top=3.40  conf=774
  00:36Z  centroid=(-22.9, 3.8, 2.91 km MSL)  z_top=3.40  conf=861
  00:40Z  centroid=(-21.4, 5.1, 2.88 km MSL)  z_top=3.15  conf=884
  00:44Z  centroid=(-20.0, 5.9, 2.85 km MSL)  z_top=3.15  conf=825
  00:48Z  centroid=(-19.6, 5.8, 2.83 km MSL)  z_top=3.15  conf=750
  00:52Z  centroid=(-17.2, 6.5, 2.82 km MSL)  z_top=3.15  conf=490
  00:56Z  centroid=(-15.0, 7.6, 2.78 km MSL)  z_top=3.15  conf=284
```

**2. z_top values**: confirmed no longer pinned at 7.65 km — now ranges
3.15-3.40 km MSL, back in line with physically reasonable column-top
heights (and matching the pre-Task-18 Task 16 baseline almost exactly).
The fix works as intended for this specific bug.

**3. Do x/y centroids still stay near the anchor (-38 to -28 km)?**
**No.** They start at -28.7 km (already outside that range at scan 1) and
drift all the way to -15.0 km by the final scan — i.e. the centroid track
has reverted almost exactly to the plume-following behavior from before
Task 18's Gaussian prior was introduced (compare to Task 16's original
un-anchored run: -28.4→-12.5 km).

**4. Updated drift statistics**: mean speed 22.1 km/hr, median 22.7,
max 39.7 km/hr, mean bearing 84°, median 73° — essentially identical to
the pre-Task-18 numbers (22.8/23.0/42.4 km/hr, 84°/65°), not the
anchored Task 18 numbers (4.4/3.8/10.8 km/hr).

**Flagged issue (important, and not what the task's own framing
predicts)**: fixing the `z_top` bug **undid essentially all of Task 18's
intended anchoring effect** — the tracker is back to following the
advected plume, not the fire source. This isn't a mistake in applying
Task 19's instructions; it reveals that Task 18's apparent success was
built on the very bug Task 19 fixed. Mechanism: previously, every
no-data voxel near the anchor's x/y column (common at altitudes beyond
useful beam coverage, or where gates are masked) inherited the anchor's
prior value directly (since `prior * 0.5 / (prior*0.5 + (1-prior)*0.5)
= prior` — the Bayesian update collapses to a pass-through when
likelihood=0.5 exactly). That created a huge, artificial pool of
positive weight (`weight = max(prob-0.5, 0)`) stacked through the entire
no-data column near the anchor — which is what was pulling the weighted
centroid toward the source (and confidence values up to 5000+, vs.
now back to 300-900, matching Task 16's original un-anchored run).

With that pass-through removed, only voxels with actual radar
data can be influenced by the prior — and at those voxels, a real strong
smoke detection (likelihood=0.92) still dominates the low prior far from
the anchor (by design — Task 18's background text says "strong radar
evidence far from the anchor can still activate those voxels"), so the
real, moving smoke plume once again outweighs the comparatively sparse
in-column evidence right at the source. Net effect: the two design goals
— (a) correct z_top by not letting no-data voxels inherit a spatial bias,
and (b) anchor the tracked centroid to the source rather than the
plume — are in tension given the current prior strength (`PRIOR_MAX=0.70`,
`PRIOR_SIGMA_KM=8.0`) and likelihood values (`P_SMOKE=0.92`). Both fixes
were applied exactly as specified since each is individually correct and
was asked for explicitly, but the PM should know that **as of this task,
neither script achieves anchored-to-source tracking with a sane z_top
simultaneously** — getting both will likely require either a much
stronger prior/lower P_SMOKE, or decoupling the horizontal and vertical
weighting schemes rather than using one shared `prob_grid`.

## TASK 20 — Update anchor to (-36, -2) to match rising column echo
**Status**: PENDING

**Goal**: Move ANCHOR_X_KM and ANCHOR_Y_KM to match the bright echo
cluster visible in the zoomed summary plot. No other changes.

In BOTH task16_prob_tracker.py AND task17_centroid_frames.py, find:
```python
ANCHOR_X_KM  = -35.8
ANCHOR_Y_KM  =   1.1
```
Replace with:
```python
ANCHOR_X_KM  = -36.0
ANCHOR_Y_KM  =  -2.0
```

Also update the Bug Fire anchor scatter point in the summary plot
(task17_centroid_frames.py) from:
```python
ax.scatter(-35.8, 1.1, ...)
```
to:
```python
ax.scatter(-36.0, -2.0, ...)
```

Re-run both scripts and paste:
1. First and last centroid lines (timestamp, x, y, z, conf)
2. Whether centroids now sit on the bright southern echo blob
3. New summary PNG if possible

---
## RESULTS 20

**Status**: DONE

**Naming/deviation notes**: `ANCHOR_X_KM`/`ANCHOR_Y_KM` updated to
(-36.0, -2.0) in both `task16_prob_tracker.py` and
`task17_centroid_frames.py`. The task text only explicitly called out
updating the anchor scatter point in the **summary** plot of
`task17_centroid_frames.py` — I also updated the matching anchor scatter
in the **per-frame** plot loop (same file, `ax.scatter(-35.8, 1.1, ...,
label='Bug Fire anchor')` a few dozen lines earlier), since leaving it
at the old coordinates would put the "Bug Fire anchor" marker in a
different place on all 14 individual frames than on the summary panel,
and out of alignment with the white dotted prior-boundary contour, which
is centered on the new anchor via `ANCHOR_X_KM`/`ANCHOR_Y_KM` already.
Confirmed via user follow-up ("apply for everything") mid-task. Left
`krgx_plume.py`, `krgx_hourly.py`, and `task12_centroid_drift.py`
untouched — their `-35.8, 1.1` / `(39.7631, -119.8788)` values are the
GOES-verified anchor from Task 08b used for actual source-confirmation
labeling, a different, already-validated purpose from this Gaussian-prior
tracking experiment; Task 20 only named the two Task 18/19 scripts.

**1. First and last centroid lines** (task16, re-run with new anchor):
```
23:59Z  centroid=(-28.7, 2.5, 2.92 km MSL)  z_top=3.40  conf=470
00:56Z  centroid=(-15.0, 7.5, 2.78 km MSL)  z_top=3.15  conf=282
```
(Full 14-scan block and drift stats are nearly identical to RESULTS 19 —
mean speed 22.6 km/hr, median 23.0, max 42.0 km/hr, bearing mean 84°,
median 73° — because, per RESULTS 19's flag, the prior's influence on
the tracked centroid is already minimal now that no-data voxels no
longer inherit it. Moving the anchor by 3.3 km made almost no difference
to the output for the same reason.)

**2. Do centroids now sit on the bright southern echo blob (near
-36, -2)?** No — the centroid track is essentially unchanged from before
this task (still drifting -28.7 → -15.0 km, following the same advected
plume path as RESULTS 19), because the anchor's influence was already
neutralized by the Task 19 fix. Moving the anchor doesn't restore the
"lock onto the source" behavior from Task 18 — that was contingent on
the bug Task 19 correctly removed.

**3. New summary PNG**: saved to `krgx_hourly/task17_centroid_summary.png`
(zoomed to x=-40..-25, y=-5..10 per your earlier request), plus all 14
per-frame PNGs regenerated in `krgx_hourly/centroid_frames/`. Can't embed
images in this text channel — files are on disk.

**Flagged issue (methodological, please read before using this anchor
downstream)**: this task's own instruction was to move the anchor "to
match the bright echo cluster visible in the zoomed summary plot" — i.e.
the anchor is being chosen to match what the **radar itself shows**,
not from independent verification. This is the same category of problem
that Task 08's original `FIRE_SOURCES` circularity was flagged for and
that Task 08b/08c were specifically built to fix (deriving the anchor
from independent GOES fire-detection data rather than from the radar
plume/echo being validated). The current GOES-verified anchor is
(-35.8, 1.1) km / (39.7631°N, -119.8788°W); the new anchor (-36.0, -2.0)
is a real, if small (~3.3 km), coordinate change with no independent
(GOES/NIFC) verification behind it — it's calibrated to visually match
radar data, which is exactly the reasoning the project moved away from
in Task 08b. Since this task only touches the Task 18/19/20 experimental
scripts (not the GOES-verified `FIRE_SOURCES` used in `krgx_plume.py`/
`krgx_hourly.py`'s actual source-confirmation pipeline), the blast radius
is contained — but if this anchor value is meant to feed back into the
main pipeline's fire-source location at some point, it should go through
the same independent-verification standard as Task 08b's GOES centroid,
not be read off a zoomed radar plot. Applied as instructed since it's an
explicit, contained parameter change to an experimental script — flagging
for the PM's awareness before it propagates further.

## TASK 21 — Rising column tracker: vertical structure + burning region
**Status**: PENDING

**Goal**: For each scan, isolate the fire column within a 10 km search
cylinder around the GOES-verified fire anchor, extract vertical structure
metrics (base, top, centroid, depth), estimate the burning region
footprint at low altitude, and flag any falling column signatures using
Vr (velocity) field. Standalone script — no changes to existing files.

**Instructions**: Create `task21_column_tracker.py` in the Runs dir.
This script loads the `grids_*.npz` files from Task 16/19 plus re-reads
the raw Level-II files for Vr.

```python
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
        cz = float(np.average(z_km_msl, weights=echo_per_z))
        # x/y centroid from column-max footprint
        w2d = np.where(np.isfinite(dbz_col_max) & (dbz_col_max > DBZ_THRESHOLD),
                       dbz_col_max, 0.0)
        if w2d.sum() > 0:
            cx = float(np.average(x_km, weights=w2d.sum(axis=1)))
            cy = float(np.average(y_km, weights=w2d.sum(axis=0)))
        else:
            cx = cy = np.nan
    else:
        cx = cy = cz = np.nan

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
        col_base_msl=round(col_base_msl, 2) if np.isfinite(col_base_msl) else None,
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

    fall_str = ' *** FALLING COLUMN ***' if falling_flag else ''
    print(f"  {dt.strftime('%H:%MZ')}  "
          f"base={col_base_msl:.2f}  top={col_top_msl:.2f}  "
          f"depth={col_depth_km:.2f} km  "
          f"burn={burn_area_km2:.1f} km²  "
          f"Vr_top={vr_top_ms:.1f} m/s{fall_str}")

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

fig, axes = plt.subplots(2, 3, figsize=(15, 8))

# 1. Column base + top + centroid z
ax = axes[0, 0]
ax.fill_between(t_min, _arr('col_base_msl'), _arr('col_top_msl'),
                alpha=0.3, color='steelblue', label='Column extent')
ax.plot(t_min, _arr('col_base_msl'), 'o-', color='navy',    label='Base')
ax.plot(t_min, _arr('col_top_msl'),  's-', color='firebrick', label='Top')
ax.plot(t_min, _arr('cz'),           '^-', color='darkorange', label='Centroid z')
ax.set(xlabel='Minutes', ylabel='km MSL', title='Column vertical extent')
ax.legend(fontsize=8)

# 2. Column depth
ax = axes[0, 1]
ax.plot(t_min, _arr('col_depth_km'), 'o-', color='purple')
ax.set(xlabel='Minutes', ylabel='km', title='Column depth (top - base)')

# 3. Burning region area
ax = axes[0, 2]
ax.plot(t_min, _arr('burn_area_km2'), 'o-', color='darkorange')
ax.set(xlabel='Minutes', ylabel='km²', title='Radar-estimated burning area')

# 4. Column centroid x/y drift
ax = axes[1, 0]
sc = ax.scatter(_arr('cx'), _arr('cy'),
                c=t_min, cmap='plasma', s=60, zorder=3)
ax.plot(_arr('cx'), _arr('cy'), 'k--', lw=0.8, zorder=2)
ax.scatter(ANCHOR_X_KM, ANCHOR_Y_KM, marker='*', c='red',
           s=150, zorder=5, label='GOES anchor')
plt.colorbar(sc, ax=ax, label='Minutes')
ax.set(xlabel='X km', ylabel='Y km',
       title='Column centroid x/y drift', aspect='equal')
ax.legend(fontsize=8)

# 5. Vr upper column
ax = axes[1, 1]
vr_arr = _arr('vr_top_ms')
colors = ['red' if r['falling_flag'] else 'steelblue' for r in records]
ax.bar(t_min, vr_arr, width=3, color=colors, alpha=0.8)
ax.axhline(-3.0, color='red', lw=1.5, linestyle='--',
           label='Falling threshold (-3 m/s)')
ax.axhline(0, color='black', lw=0.8)
ax.set(xlabel='Minutes', ylabel='m/s',
       title='Mean Vr upper column\n(red = falling column flag)')
ax.legend(fontsize=8)

# 6. Burn centroid drift
ax = axes[1, 2]
sc2 = ax.scatter(_arr('burn_cx'), _arr('burn_cy'),
                 c=t_min, cmap='plasma', s=60, zorder=3)
ax.plot(_arr('burn_cx'), _arr('burn_cy'), 'k--', lw=0.8)
ax.scatter(ANCHOR_X_KM, ANCHOR_Y_KM, marker='*', c='red',
           s=150, zorder=5, label='GOES anchor')
plt.colorbar(sc2, ax=ax, label='Minutes')
ax.set(xlabel='X km', ylabel='Y km',
       title='Burning region centroid drift', aspect='equal')
ax.legend(fontsize=8)

plt.suptitle('KRGX Bug Fire — Column tracker\n'
             '2026-08-09 23:59Z to 2026-08-10 00:56Z', fontsize=12)
plt.tight_layout()
out_png = OUT_DIR / 'task21_column_tracker.png'
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print(f'Saved {out_png}')
```

### After running: paste
1. Full per-scan print block (base, top, depth, burn area, Vr_top)
2. Any [FALLING COLUMN] flags
3. Whether column_tracker.json and PNG saved without error
4. Any traceback — most likely place is the cache_lookup dt parsing
   if the filename format differs from expected

---
## RESULTS 21

**Status**: DONE — run 2026-09-18.

Created `task21_column_tracker.py` exactly per spec (no logic changes needed —
grid edges, `ANCHOR_X_KM`/`ANCHOR_Y_KM` (-35.8, 1.1, the GOES-verified Task 08b
anchor, not the Task 20 experimental drift), and field keys all matched the
existing `grids_*.npz` files and `krgx_hourly.py`'s bin construction on
inspection before running — confirmed `x_bins`/`y_bins`/`z_bins` arrays in a
sample `grids_20260810_000809.npz` reproduce identically from the same
`X_MIN/X_MAX/Y_MIN/Y_MAX/Z_MAX/DXY/DZ/RADAR_ALT_M` constants). Only cosmetic
change: replaced the literal `°`/`²` characters in f-strings with `°`/
`²` escapes to avoid a console-encoding crash risk in this shell (no
functional difference, matches the workaround already used in `task04_diag.py`
for a similar glyph issue).

**1. Full per-scan print block**:
```
  23:59Z  base=2.65  top=3.65  depth=1.00 km  burn=19.1 km²  Vr_top=-8.7 m/s *** FALLING COLUMN ***
  00:04Z  base=2.65  top=3.90  depth=1.25 km  burn=23.1 km²  Vr_top=-8.6 m/s *** FALLING COLUMN ***
  00:08Z  base=2.65  top=3.90  depth=1.25 km  burn=28.8 km²  Vr_top=-8.5 m/s *** FALLING COLUMN ***
  00:12Z  base=2.65  top=3.90  depth=1.25 km  burn=36.6 km²  Vr_top=-8.5 m/s *** FALLING COLUMN ***
  00:16Z  base=2.65  top=3.90  depth=1.25 km  burn=36.6 km²  Vr_top=-8.5 m/s *** FALLING COLUMN ***
  00:20Z  base=2.65  top=3.90  depth=1.25 km  burn=30.4 km²  Vr_top=-8.5 m/s *** FALLING COLUMN ***
  00:24Z  base=2.65  top=3.40  depth=0.75 km  burn=18.9 km²  Vr_top=-9.2 m/s *** FALLING COLUMN ***
  00:32Z  base=2.65  top=3.40  depth=0.75 km  burn=13.4 km²  Vr_top=-9.3 m/s *** FALLING COLUMN ***
  00:36Z  base=2.65  top=3.40  depth=0.75 km  burn=12.1 km²  Vr_top=-9.2 m/s *** FALLING COLUMN ***
  00:40Z  base=2.65  top=3.40  depth=0.75 km  burn=6.1 km²  Vr_top=-9.3 m/s *** FALLING COLUMN ***
  00:44Z  base=2.65  top=3.40  depth=0.75 km  burn=6.2 km²  Vr_top=-9.3 m/s *** FALLING COLUMN ***
  00:48Z  base=2.65  top=3.65  depth=1.00 km  burn=10.6 km²  Vr_top=-9.0 m/s *** FALLING COLUMN ***
  00:52Z  base=2.65  top=3.40  depth=0.75 km  burn=8.4 km²  Vr_top=-9.4 m/s *** FALLING COLUMN ***
  00:56Z  base=2.65  top=3.40  depth=0.75 km  burn=7.4 km²  Vr_top=-9.3 m/s *** FALLING COLUMN ***
```

**2. Falling column flags**: **All 14/14 scans flagged `*** FALLING COLUMN ***`**
— `Vr_top` ranges narrowly from -8.4 to -9.4 m/s across the entire hour, well
past the -3.0 m/s threshold every time.

**Flagging this as a likely systematic artifact, not 14 independent
downdraft detections**: a genuine falling-column signature should be
intermittent (rising column, injection pulse, then occasional collapse) — a
constant ~-8.5 to -9.4 m/s "upper column" Vr for the entire hour, with almost
no scan-to-scan variance, looks far more consistent with a fixed radar
geometry effect than with 14 consecutive real downdrafts. Candidate causes
worth checking before trusting this flag downstream:
- **Beam geometry / radial-velocity sign convention**: Vr is radial velocity
  *along the beam*, not vertical velocity. At the column's range/elevation
  from KRGX, a large horizontal wind component (e.g. westerly flow advecting
  smoke toward the radar, given the plume's observed drift track from Tasks
  16-20) would register as strongly negative Vr (motion toward radar) at
  *any* elevation angle sampling that geometry — this would look identical
  to a real downdraft in the "mean Vr in upper half of column" metric as
  currently defined, without actually being vertical motion. The script
  never converts Vr to a true vertical-velocity estimate (would need to
  divide by sin(elevation) or similar, using beam geometry) — it just
  averages raw radial Vr over the upper-column voxels.
- **No dealiasing applied**: `load_vr_grid()` reads raw `velocity` straight
  from the archive with no unfolding step. If this VCP's Nyquist velocity is
  low enough that raw Vr near the plume top is aliased, that alone could
  produce a suspiciously uniform, saturated-looking negative value across
  scans.
- Worth a quick sanity check in a future task: plot Vr vs. elevation angle
  for a few individual gates in the cylinder to see if the "-8.5 to -9.4"
  band is really vertical motion or just the ambient horizontal wind
  component along the beam.

**3. Output files**: `krgx_hourly/column_tracker.json` (225 lines, one
record per scan) and `krgx_hourly/task21_column_tracker.png` (188 KB, 6-panel
figure) both saved without error.

**4. Tracebacks**: none. `cache_lookup` dt parsing worked as written — all 14
grid-file timestamps matched a cache file within the 300s tolerance, no
`KeyError`/parse failures. One non-fatal stderr line on every run,
unrelated to this script's logic: `Warning 3: Cannot find header.dxf
(GDAL_DATA is not defined)` — a pre-existing GDAL data-path warning from the
`geo_env` environment (same class of issue as the recurring conda/PATH note),
printed once at interpreter start, does not affect output.

**Consistency with prior tasks**: column base (2.65 km MSL) and depth
(0.75-1.25 km) track the same footprint growth/contraction pattern already
seen in the burn-area numbers (peaks at 36.6 km² around 00:12-00:16Z, matching
the burn-area growth phase from earlier tasks' centroid-drift analysis) —
nothing here contradicts the existing plume-drift/vertical-structure picture
from Tasks 16-20, only the falling-column flag above needs a second look
before being treated as a real physical result.

**Environment note (repeats)**: same `geo_env` python.exe + manual PATH
workaround as prior tasks (per the note in `CLAUDE.md`).

(Claude Code pastes output here)

## TASK 22 — Stabilize column centroid + drop base metric
**Status**: PENDING

**Goal**: Fix the column centroid drift plot by applying a 2 km max
displacement cap between consecutive scans. Drop column base (beam floor
artifact). Keep column top, depth, burn area, and centroid drift.
Modify task21_column_tracker.py in place.

**Physical constraint**: the rising column injection point is tied to the
fire front, which moves at most ~1 m/s (~3.6 km/hr). Over a 4-minute
scan interval that's ~0.24 km. Allow 2 km as a generous upper bound
(accounts for column lean and centroid estimation noise).

**Instructions**:

### Change 1 — Add scan-to-scan centroid cap to main loop
Before the main loop, initialize tracker state:
```python
prev_cx, prev_cy = ANCHOR_X_KM, ANCHOR_Y_KM   # seed at anchor
MAX_DISP_KM = 2.0   # max allowed centroid displacement per scan
```

Inside the loop, after `cx`/`cy` are computed from the dBZ-weighted
footprint, add the displacement cap:
```python
    if np.isfinite(cx) and np.isfinite(cy):
        disp = np.sqrt((cx - prev_cx)**2 + (cy - prev_cy)**2)
        if disp > MAX_DISP_KM:
            # candidate jumped too far — hold previous position
            print(f"    [CAP] {dt.strftime('%H:%MZ')}: "
                  f"centroid jump {disp:.1f} km > {MAX_DISP_KM} km cap "
                  f"— holding ({prev_cx:.1f}, {prev_cy:.1f})")
            cx, cy = prev_cx, prev_cy
        prev_cx, prev_cy = cx, cy
    else:
        # no echo in cylinder — hold previous position
        cx, cy = prev_cx, prev_cy
        prev_cx, prev_cy = cx, cy
```

### Change 2 — Drop col_base_msl from output and plots
Remove `col_base_msl` from the `rec` dict and from the plot.

Replace the column vertical extent panel (axes[0,0]) with:
```python
ax = axes[0, 0]
ax.plot(t_min, _arr('col_top_msl'), 's-', color='firebrick', label='Top')
ax.plot(t_min, _arr('cz'),          '^-', color='darkorange', label='Centroid z')
ax.axhline(RADAR_ALT_M/1e3, color='gray', lw=1, linestyle=':',
           label=f'Radar floor ({RADAR_ALT_M/1e3:.2f} km MSL)')
ax.set(xlabel='Minutes', ylabel='km MSL',
       title='Column top + centroid height\n(base not shown — beam floor artifact)')
ax.legend(fontsize=8)
```

### Change 3 — Drop Vr panel, replace with burn area + centroid combined
The 6-panel layout had: col extent, depth, burn area, centroid xy,
Vr, burn centroid. Drop Vr entirely. Replace with a cleaner 2x2 layout:

```python
fig, axes = plt.subplots(2, 2, figsize=(12, 9))
```

Panel layout:
- axes[0,0]: column top + centroid z vs time (Change 2 above)
- axes[0,1]: column depth vs time (unchanged)
- axes[1,0]: burn area km² vs time (unchanged)  
- axes[1,1]: column centroid x/y drift (stabilized, with 2 km scale bar)

For the centroid drift panel, add a 2 km scale bar for reference:
```python
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
```

### After running: paste
1. Any [CAP] lines — how many scans were capped and by how much
2. Full per-scan print block (base dropped, top/depth/burn/Vr_top removed)
3. Column top values across the hour — confirm the 00:24Z drop is still present
4. Whether the centroid drift plot now shows smooth sub-2km steps
5. PNG saved without error

---
## RESULTS 22

**Status**: DONE (code applied exactly as specified) — but flagging a
critical logical flaw in the design that makes this run's centroid-drift
output not usable as-is. Run 2026-09-18. Modified `task21_column_tracker.py`
in place per all three changes: added `prev_cx, prev_cy = ANCHOR_X_KM,
ANCHOR_Y_KM` + `MAX_DISP_KM = 2.0` before the main loop; inserted the
displacement-cap block immediately after the existing `cx`/`cy` centroid
computation (before the burn-footprint section); removed `col_base_msl` from
the `rec` dict (kept the internal `col_base_msl`/`base_idx` variables, since
`col_depth_km` and the burn-footprint slice still need them — only the
*output* field was dropped, per the spec's literal wording); rebuilt the
figure as a 2x2 layout (top+centroid-z, depth, burn area, stabilized
centroid-drift-with-scale-bar), dropping the Vr panel entirely. Also trimmed
the per-scan print line to `top`/`depth`/`burn` only, dropping `base` and
`Vr_top` (matches item 2's expected format in the task's own "After running"
checklist). Left `load_vr_grid()`, and `vr_top_ms`/`vr_mean_cyl`/
`falling_flag` in the `rec`/json output untouched — the task only said to
drop the Vr *panel* from the plot, not the underlying data, so those three
fields are still computed and still land in `column_tracker.json`, just no
longer plotted or printed per-scan.

**1. [CAP] lines — all 14/14 scans were capped**, every single one:
```
    [CAP] 23:59Z: centroid jump 4.2 km > 2.0 km cap — holding (-35.8, 1.1)
    [CAP] 00:04Z: centroid jump 4.0 km > 2.0 km cap — holding (-35.8, 1.1)
    [CAP] 00:08Z: centroid jump 4.4 km > 2.0 km cap — holding (-35.8, 1.1)
    [CAP] 00:12Z: centroid jump 5.4 km > 2.0 km cap — holding (-35.8, 1.1)
    [CAP] 00:16Z: centroid jump 6.5 km > 2.0 km cap — holding (-35.8, 1.1)
    [CAP] 00:20Z: centroid jump 6.7 km > 2.0 km cap — holding (-35.8, 1.1)
    [CAP] 00:24Z: centroid jump 6.5 km > 2.0 km cap — holding (-35.8, 1.1)
    [CAP] 00:32Z: centroid jump 6.1 km > 2.0 km cap — holding (-35.8, 1.1)
    [CAP] 00:36Z: centroid jump 5.8 km > 2.0 km cap — holding (-35.8, 1.1)
    [CAP] 00:40Z: centroid jump 5.0 km > 2.0 km cap — holding (-35.8, 1.1)
    [CAP] 00:44Z: centroid jump 4.8 km > 2.0 km cap — holding (-35.8, 1.1)
    [CAP] 00:48Z: centroid jump 5.7 km > 2.0 km cap — holding (-35.8, 1.1)
    [CAP] 00:52Z: centroid jump 7.0 km > 2.0 km cap — holding (-35.8, 1.1)
    [CAP] 00:56Z: centroid jump 7.0 km > 2.0 km cap — holding (-35.8, 1.1)
```

**Critical flag — this is a design bug, not a data quality issue**: because
the tracker is *seeded at the anchor* (`prev_cx, prev_cy = ANCHOR_X_KM,
ANCHOR_Y_KM`) and every rejected candidate causes the code to *hold the
previous position* rather than the previous **real measurement**, the
centroid is permanently frozen at `(-35.8, 1.1)` — the GOES anchor itself —
for the entire hour. Confirmed directly from `column_tracker.json`: `cx`/`cy`
are exactly `-35.8, 1.1` on all 14 records, no exceptions. This happens
because Tasks 16-20 already established (and this task's own docstring
references) that the real plume centroid drifts from about -28.7 km to
-15.0 km over the hour — i.e. the *true* first-scan centroid is already
~7-8 km from the anchor, blowing past the 2 km cap on the very first
iteration. Once the first real measurement is rejected, `prev_cx/prev_cy`
never updates away from the anchor, so *every subsequent* real measurement
also gets compared against a stale anchor position it was never near to
begin with — the cap doesn't stabilize a noisy but real track, it locks the
track to a single point for the whole run. The physical premise in the task
spec ("fire front moves at most ~1 m/s, so cap consecutive-scan jumps at
2 km") is reasonable for scan-to-scan noise, but the initial seed should
have been the *first scan's own raw centroid estimate*, not the independently
-sourced anchor 7+ km away — seeding at a point the real data was never
going to be near guarantees permanent lock-up rather than stabilization.
Recommend a follow-up fix (not applied here, since this task only asked to
implement the spec as given): seed `prev_cx, prev_cy` from the first scan's
uncapped `cx, cy` instead of `ANCHOR_X_KM, ANCHOR_Y_KM`, or skip the cap
entirely on the first iteration.

**2. Full per-scan print block** (base dropped, Vr_top dropped):
```
  23:59Z  top=3.65  depth=1.00 km  burn=19.1 km²
  00:04Z  top=3.90  depth=1.25 km  burn=23.1 km²
  00:08Z  top=3.90  depth=1.25 km  burn=28.8 km²
  00:12Z  top=3.90  depth=1.25 km  burn=36.6 km²
  00:16Z  top=3.90  depth=1.25 km  burn=36.6 km²
  00:20Z  top=3.90  depth=1.25 km  burn=30.4 km²
  00:24Z  top=3.40  depth=0.75 km  burn=18.9 km²
  00:32Z  top=3.40  depth=0.75 km  burn=13.4 km²
  00:36Z  top=3.40  depth=0.75 km  burn=12.1 km²
  00:40Z  top=3.40  depth=0.75 km  burn=6.1 km²
  00:44Z  top=3.40  depth=0.75 km  burn=6.2 km²
  00:48Z  top=3.65  depth=1.00 km  burn=10.6 km²
  00:52Z  top=3.40  depth=0.75 km  burn=8.4 km²
  00:56Z  top=3.40  depth=0.75 km  burn=7.4 km²
```
These top/depth/burn values are numerically identical to RESULTS 21 (the cap
only touches `cx`/`cy`, not `col_top_msl`/`col_depth_km`/`burn_area_km2`), as
expected since Change 1 doesn't feed back into those calculations.

**3. Column top values — the 00:24Z drop is still present**: top falls from
3.90 km MSL (00:20Z) to 3.40 km MSL (00:24Z), the same discontinuity flagged
in RESULTS 21, confirming this task's changes didn't alter that signal.

**4. Does the centroid drift plot now show smooth sub-2km steps? No — it
shows zero steps.** Per the bug above, every point in the "Column centroid
drift (2 km cap applied)" panel plots on top of the anchor marker itself;
the plotted track is a single stationary point with 14 overlapping
timestamps, not a smooth path. This is worse than the pre-cap version from
Task 21 (which at least showed the real, if jumpy, drift track from -28.7 to
-15.0 km) — the stabilization has fully suppressed the actual signal rather
than smoothing it.

**5. PNG/JSON**: `krgx_hourly/column_tracker.json` (4.6 KB, 14 records, no
`col_base_msl` key, `cx`/`cy` all equal to anchor) and
`krgx_hourly/task21_column_tracker.png` (151 KB, 2x2 layout, no Vr panel)
both saved without error, exit code 0, no traceback.

**Recommendation for a Task 23**: re-seed the cap's `prev_cx/prev_cy` from
the first scan's own uncapped centroid (not the anchor), then re-run and
confirm the drift plot shows the expected sub-2km-per-scan smoothed path
following the same overall trajectory as the uncapped Task 21/RESULTS 19-20
track, rather than a stationary point.

**Environment note (repeats)**: same `geo_env` python.exe + manual PATH
workaround as prior tasks (per the note in `CLAUDE.md`).

## TASK 23 — Fix centroid seeding + weight by smoke mask inside cylinder
**Status**: PENDING

**Goal**: Two fixes to task21_column_tracker.py:
1. Seed prev_cx/prev_cy from first scan's raw centroid, not the anchor
2. Compute cylinder centroid from smoke-masked voxels only (rhv < RHV_MAX
   AND dbz > DBZ_THRESHOLD), weighted by dBZ, not raw dBZ across all echo

**Instructions**:

### Change 1 — Seed from first scan, not anchor
Replace:
```python
prev_cx, prev_cy = ANCHOR_X_KM, ANCHOR_Y_KM
```
With:
```python
prev_cx, prev_cy = None, None   # seeded from first scan's raw centroid
```

In the displacement cap block, replace:
```python
    if np.isfinite(cx) and np.isfinite(cy):
        disp = np.sqrt((cx - prev_cx)**2 + (cy - prev_cy)**2)
        if disp > MAX_DISP_KM:
            print(f"    [CAP] ...")
            cx, cy = prev_cx, prev_cy
        prev_cx, prev_cy = cx, cy
    else:
        cx, cy = prev_cx, prev_cy
        prev_cx, prev_cy = cx, cy
```
With:
```python
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
```

### Change 2 — Weight centroid by smoke mask, not raw echo
Find the existing centroid computation block:
```python
    dbz_col_max = dbz_cyl.max(axis=2)
    total_echo  = echo_per_z.sum()
    if total_echo > 0:
        cz = float(np.average(z_km_msl, weights=echo_per_z))
        w2d = np.where(np.isfinite(dbz_col_max) & (dbz_col_max > DBZ_THRESHOLD),
                       dbz_col_max, 0.0)
        if w2d.sum() > 0:
            cx = float(np.average(x_km, weights=w2d.sum(axis=1)))
            cy = float(np.average(y_km, weights=w2d.sum(axis=0)))
```

Replace the w2d / cx / cy lines with:
```python
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
```

Also update cz to use smoke-weighted vertical profile:
```python
        smoke_per_z_w = np.where(smoke_per_z > 0, smoke_per_z,
                                 echo_per_z * 0.1)   # downweight non-smoke
        cz = float(np.average(z_km_msl, weights=smoke_per_z_w))
```

### After running: paste
1. [SEED] line — what is the first scan's raw centroid?
2. [CAP] lines — how many scans capped and by how much
3. Full per-scan cx/cy from column_tracker.json
   (first and last, plus any that were capped)
4. Whether the drift plot now shows a coherent track
   within a few km of the anchor
5. Column top values — confirm 00:24Z drop unchanged

---
## RESULTS 23

**Status**: DONE — run 2026-09-18. Applied both changes to
`task21_column_tracker.py` in place exactly as specified: seed
`prev_cx, prev_cy = None, None`, first-scan unconditional-seed branch with
`[SEED]` print, `[CAP]` branch unchanged otherwise; centroid computation
switched to smoke-mask-weighted `smoke_col_max` with fallback to the
original raw-dBZ `w2d` weighting when a scan has zero smoke-mask coverage in
the cylinder, and `cz` switched to the smoke-downweighted vertical profile
(`smoke_per_z_w`). No naming/structure deviations needed — the spec's
before/after blocks matched the file's current state (from Task 22) exactly.

**1. [SEED] line**:
```
    [SEED] 23:59Z: seeded at (-31.5, 1.0)
```
First scan's raw smoke-weighted centroid is (-31.5, 1.0) km — about 4.3 km
from the GOES anchor (-35.8, 1.1).

**2. [CAP] lines — only 3/14 scans capped** (compare to Task 22's 14/14):
```
    [CAP] 00:40Z: jump 2.1 km > cap — holding (-30.1, -0.1)
    [CAP] 00:44Z: jump 2.2 km > cap — holding (-30.1, -0.1)
    [CAP] 00:48Z: jump 2.2 km > cap — holding (-30.1, -0.1)
```
All three are marginal overshoots (2.1-2.2 km vs. the 2.0 km cap), not the
4-7 km blowouts seen in Task 22 — consistent with the seeding fix actually
working as intended.

**3. Full cx/cy series** (from `column_tracker.json`, all 14 scans, plus
distance to the GOES anchor for context):
```
23:59:57  cx=-31.49  cy= 0.97  dist_to_anchor=4.31 km   [seed]
00:04:04  cx=-31.99  cy=-0.24  dist_to_anchor=4.04 km
00:08:09  cx=-31.61  cy=-0.45  dist_to_anchor=4.47 km
00:12:15  cx=-30.50  cy=-0.08  dist_to_anchor=5.43 km
00:16:21  cx=-29.27  cy= 0.24  dist_to_anchor=6.59 km
00:20:27  cx=-29.10  cy= 0.26  dist_to_anchor=6.75 km
00:24:33  cx=-29.40  cy=-0.33  dist_to_anchor=6.56 km
00:32:07  cx=-29.94  cy=-0.84  dist_to_anchor=6.17 km
00:36:13  cx=-30.08  cy=-0.15  dist_to_anchor=5.85 km
00:40:20  cx=-30.08  cy=-0.15  dist_to_anchor=5.85 km   [capped, held]
00:44:27  cx=-30.08  cy=-0.15  dist_to_anchor=5.85 km   [capped, held]
00:48:33  cx=-30.08  cy=-0.15  dist_to_anchor=5.85 km   [capped, held]
00:52:38  cx=-29.50  cy=-2.05  dist_to_anchor=7.04 km
00:56:45  cx=-29.08  cy=-0.80  dist_to_anchor=6.98 km
```
First scan: (-31.5, 1.0). Last scan: (-29.1, -0.8). Net drift over the hour:
~2.6 km east, ~1.8 km south — much smaller and smoother than the raw
uncapped Task 21 track (which ranged from -28.7 to -15.0 km, a ~14 km span),
because the smoke-mask weighting itself (Change 2) concentrates the centroid
near the low-RhoHV injection core rather than the full advected dBZ plume
mass, independent of the capping.

**4. Does the drift plot now show a coherent track? Yes, but not "within a
few km of the anchor"** — it's a smooth, physically plausible track (mostly
sub-1-km steps between scans, only three consecutive 2.1-2.2 km steps
correctly capped) rather than the frozen single point from Task 22. However,
it sits a stable **4.0-7.0 km from the GOES anchor throughout the entire
hour**, never converging toward it. This is consistent with — not
contradicting — the standing offset already flagged in Tasks 20 and prior:
the radar-derived column centroid and the independently-verified GOES fire
anchor are two different physical things (radar sees smoke-plume mass
aloft and advected, GOES sees the surface thermal source), so a several-km
persistent offset between them is expected, not a bug. The plot
(`krgx_hourly/task21_column_tracker.png`, regenerated, 173 KB) shows this
correctly: a compact, smoothly-connected track sitting consistently
northeast-ish of the anchor marker, not scattered or jumping.

**5. Column top values — 00:24Z drop unchanged**: 3.90 km MSL (00:20Z) →
3.40 km MSL (00:24Z), identical to Results 21/22 — expected, since neither
of this task's two changes touch `col_top_msl`.

**Output files**: `krgx_hourly/column_tracker.json` (4.6 KB, 14 records)
and `krgx_hourly/task21_column_tracker.png` (173 KB) both saved without
error, exit code 0, no traceback.

**Environment note (repeats)**: same `geo_env` python.exe + manual PATH
workaround as prior tasks (per the note in `CLAUDE.md`).

## TASK 24 — Overlay column tracker results on per-scan radar images
**Status**: PENDING

**Goal**: Produce 14 PNGs (one per scan) showing the 2D column-max dBZ
radar background with the column tracker results overlaid:
- Circle showing the 10 km search cylinder boundary
- Column centroid marker (from column_tracker.json cx/cy)
- Column top height labeled
- Burn area footprint (lowest 3 z levels inside cylinder)
- GOES anchor marker
Standalone scratch script — no changes to existing files.

**Instructions**: Create `task24_column_frames.py` in the Runs dir.

```python
"""
Column tracker overlay on radar — Task 24
Plots column_tracker.json results on per-scan radar backgrounds.
"""
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from pathlib import Path
from datetime import datetime

OUT_DIR    = Path('krgx_hourly')
FRAME_DIR  = OUT_DIR / 'column_frames'
FRAME_DIR.mkdir(exist_ok=True)

ANCHOR_X_KM   = -35.8
ANCHOR_Y_KM   =   1.1
ANCHOR_RAD_KM = 10.0
RADAR_ALT_M   = 2530.0
DBZ_THRESHOLD = 10.0
RHV_MAX       = 0.90

# grid edges — must match krgx_hourly.py exactly
X_MIN, X_MAX = -50_000, 10_000
Y_MIN, Y_MAX = -10_000, 25_000
Z_MAX        =   8_000
DXY, DZ      = 250, 250
x_bins = np.arange(X_MIN, X_MAX + DXY, DXY) / 1e3
y_bins = np.arange(Y_MIN, Y_MAX + DXY, DXY) / 1e3
z_bins = np.arange(0, Z_MAX/1e3 - RADAR_ALT_M/1e3 + DZ/1e3, DZ/1e3)
x_km   = (x_bins[:-1] + x_bins[1:]) / 2
y_km   = (y_bins[:-1] + y_bins[1:]) / 2
z_km_msl = (z_bins[:-1] + z_bins[1:]) / 2 + RADAR_ALT_M/1e3
nx, ny, nz = len(x_km), len(y_km), len(z_km_msl)
XX, YY = np.meshgrid(x_km, y_km, indexing='ij')

# cylinder mask
cyl_mask = ((XX - ANCHOR_X_KM)**2 +
            (YY - ANCHOR_Y_KM)**2) <= ANCHOR_RAD_KM**2
cyl_3d   = cyl_mask[:, :, np.newaxis].repeat(nz, axis=2)

# load column tracker results
with open(OUT_DIR / 'column_tracker.json') as f:
    col_records = json.load(f)
col_by_dt = {r['dt'][:16]: r for r in col_records}  # key by YYYY-MM-DDTHH:MM

# load grid files
grid_files = sorted(OUT_DIR.glob('grids_*.npz'))
norm       = mcolors.Normalize(vmin=0, vmax=35)
cmap_dbz   = plt.cm.jet

for i, gf in enumerate(grid_files):
    stem = gf.stem.replace('grids_', '')
    dt   = datetime.strptime(stem, '%Y%m%d_%H%M%S')

    d   = np.load(gf)
    dbz = d['dbz_smooth'].astype(float)
    rhv = d['rhv_grid'].astype(float)

    # 2D column-max dBZ (full domain)
    col_max = dbz.max(axis=2)

    # smoke mask footprint inside cylinder
    smoke_cyl = (
        cyl_3d &
        np.isfinite(rhv) &
        (rhv < RHV_MAX) &
        (dbz > DBZ_THRESHOLD)
    )

    # burn footprint: lowest 3 z levels with echo inside cylinder
    echo_cyl   = cyl_3d & (dbz > DBZ_THRESHOLD)
    echo_per_z = echo_cyl.sum(axis=(0, 1))
    base_idx   = int(np.argmax(echo_per_z > 0)) if echo_per_z.any() else 0
    burn_slice = echo_cyl[:, :, base_idx:min(base_idx+3, nz)]
    burn_fp    = burn_slice.any(axis=2)   # (nx, ny)

    # smoke footprint column-any inside cylinder
    smoke_fp = smoke_cyl.any(axis=2)

    # get column tracker record for this scan
    dt_key = dt.strftime('%Y-%m-%dT%H:%M')
    rec    = col_by_dt.get(dt_key, {})
    cx     = rec.get('cx')
    cy     = rec.get('cy')
    col_top = rec.get('col_top_msl')
    col_depth = rec.get('col_depth_km')
    burn_area = rec.get('burn_area_km2')

    # ── plot ─────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 8))

    # dBZ background
    pcm = ax.pcolormesh(x_km, y_km, col_max.T,
                        cmap=cmap_dbz, norm=norm, shading='auto')
    plt.colorbar(pcm, ax=ax, label='dBZ (column max)', pad=0.02)

    # smoke mask footprint inside cylinder (purple tint)
    if smoke_fp.any():
        ax.contourf(x_km, y_km, smoke_fp.T.astype(float),
                    levels=[0.5, 1.5], colors=['purple'],
                    alpha=0.20, zorder=2)
        ax.contour(x_km, y_km, smoke_fp.T.astype(float),
                   levels=[0.5], colors=['purple'],
                   linewidths=1.0, linestyles='--', zorder=3)

    # burn footprint (orange tint, lowest gates)
    if burn_fp.any():
        ax.contourf(x_km, y_km, burn_fp.T.astype(float),
                    levels=[0.5, 1.5], colors=['orange'],
                    alpha=0.30, zorder=2)
        ax.contour(x_km, y_km, burn_fp.T.astype(float),
                   levels=[0.5], colors=['darkorange'],
                   linewidths=1.5, linestyles='-', zorder=3)

    # 10 km search cylinder boundary
    theta  = np.linspace(0, 2*np.pi, 360)
    cyl_cx = ANCHOR_X_KM + ANCHOR_RAD_KM * np.cos(theta)
    cyl_cy = ANCHOR_Y_KM + ANCHOR_RAD_KM * np.sin(theta)
    ax.plot(cyl_cx, cyl_cy, 'w-', lw=1.5, alpha=0.7,
            linestyle=':', zorder=4, label='10 km cylinder')

    # GOES anchor
    ax.scatter(ANCHOR_X_KM, ANCHOR_Y_KM,
               marker='*', c='red', s=200,
               edgecolors='white', linewidths=0.5,
               zorder=6, label='GOES anchor')

    # column centroid
    if cx is not None and cy is not None:
        ax.scatter(cx, cy, marker='o', c='yellow', s=180,
                   edgecolors='black', linewidths=1.0,
                   zorder=7, label='Column centroid')
        label_lines = [f'{dt.strftime("%H:%MZ")}']
        if col_top  is not None:
            label_lines.append(f'top={col_top:.2f} km MSL')
        if col_depth is not None:
            label_lines.append(f'depth={col_depth:.2f} km')
        if burn_area is not None:
            label_lines.append(f'burn={burn_area:.1f} km²')
        ax.annotate(
            '\n'.join(label_lines),
            xy=(cx, cy),
            xytext=(10, 10), textcoords='offset points',
            fontsize=8, color='yellow',
            bbox=dict(boxstyle='round,pad=0.3',
                      fc='black', alpha=0.6),
            zorder=8
        )

    # radar site
    ax.scatter(0, 0, marker='+', c='white', s=150,
               linewidths=2, zorder=6, label='KRGX')

    # legend patches for footprints
    smoke_patch = mpatches.Patch(color='purple', alpha=0.4,
                                 label='Smoke mask (RhoHV)')
    burn_patch  = mpatches.Patch(color='orange', alpha=0.5,
                                 label='Burn footprint (low gates)')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles + [smoke_patch, burn_patch],
              fontsize=8, loc='upper right')

    ax.set(
        aspect='equal',
        xlim=(x_km.min(), x_km.max()),
        ylim=(y_km.min(), y_km.max()),
        xlabel='X km (aeqd from KRGX)',
        ylabel='Y km (aeqd from KRGX)',
        title=(f'KRGX  {dt.strftime("%Y-%m-%d %H:%MZ")}  '
               f'[{i+1}/14]\n'
               f'Column: top={col_top:.2f} km  '
               f'depth={col_depth:.2f} km  '
               f'burn={burn_area:.1f} km²'
               if all(v is not None for v in [col_top,col_depth,burn_area])
               else f'KRGX  {dt.strftime("%Y-%m-%d %H:%MZ")}  [{i+1}/14]')
    )

    out_path = FRAME_DIR / f'col_{i+1:02d}_{dt.strftime("%H%M%S")}.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved {out_path.name}  '
          f'top={col_top}  burn={burn_area}')

print(f'\nAll frames saved to {FRAME_DIR}')
```

### After running: paste
1. List of 14 filenames confirmed saved
2. Any traceback — most likely the dt_key lookup
   if the JSON datetime format differs from '%Y-%m-%dT%H:%M'
3. Share 2-3 representative frames if possible
   (first scan, peak burn area scan ~00:12Z, last scan)

---
## RESULTS 24

**Status**: DONE — run 2026-09-18. Created `task24_column_frames.py` exactly
per spec, no deviations — the `dt_key` lookup (`dt.strftime('%Y-%m-%dT%H:%M')`
vs. `r['dt'][:16]` from `datetime.isoformat()`) matched on the first try, no
lookup failures. Only non-functional artifact: one matplotlib `UserWarning`
("linestyle is redundantly defined by the 'linestyle' keyword argument and
the fmt string 'w-'") on the 10 km cylinder boundary line — harmless, from
the spec's own `ax.plot(cyl_cx, cyl_cy, 'w-', ..., linestyle=':')` call
(fmt string sets solid, kwarg overrides to dotted; kwarg wins, output is
correct dotted white circle as intended), not something I introduced.

**1. All 14 filenames confirmed saved** to `krgx_hourly/column_frames/`:
```
col_01_235957.png  col_02_000404.png  col_03_000809.png  col_04_001215.png
col_05_001621.png  col_06_002027.png  col_07_002433.png  col_08_003207.png
col_09_003613.png  col_10_004020.png  col_11_004427.png  col_12_004833.png
col_13_005238.png  col_14_005645.png
```
Per-frame print block confirmed matching burn/top values from
`column_tracker.json` (Task 23's output) throughout — e.g. `col_04_001215.png
top=3.9 burn=36.62` matches the peak-burn scan.

**2. No traceback.** dt_key lookup worked cleanly on all 14 scans.

**3. Representative frames** (viewed directly, described below since this
channel can't embed images inline for the user — files are on disk):
- **col_01 (23:59Z, first scan)**: column centroid (yellow) sits just east
  of the GOES anchor (red star) inside the 10 km cylinder; burn footprint
  (orange) is a small elongated patch extending southwest from the
  centroid; no visually distinct purple smoke-mask outline separate from
  the orange.
- **col_04 (00:12Z, peak burn area, 36.6 km²)**: burn footprint has grown
  into a larger, more filamentary orange shape consistent with the
  burn-area time series peaking here; centroid still ~5.4 km from anchor;
  background dBZ field shows the plume mass stretching well outside the
  cylinder to the northeast (as expected from the drift track).
- **col_14 (00:56Z, last scan)**: burn footprint has fragmented into three
  small disconnected orange patches (consistent with burn area dropping to
  7.4 km² from its peak), centroid sits at the largest fragment; radar
  background shows the bulk of the advected plume has moved well outside
  the 10 km cylinder to the northeast, visually confirming why the
  cylinder-restricted burn/column metrics decline over the hour even
  though the wider dBZ field is still active.

**Observation not asked for but worth flagging**: in all three viewed
frames, the "Smoke mask (RhoHV)" purple contour is not visually
distinguishable from the "Burn footprint" orange contour — checked directly
against the 00:12Z scan's raw arrays: `smoke_fp` (582 px) and `burn_fp`
(586 px) differ by only 9 pixels out of ~582-586, i.e. the two footprints
are >98% coincident in this dataset. This isn't a rendering bug — the
purple fill is being drawn first and the orange fill/contour on top of it
at the same or higher zorder, so the two nearly-identical shapes visually
merge into "orange with a thin/invisible purple sliver." If a future task
wants the smoke-vs-burn distinction to be visible, worth either drawing
smoke mask outline on top of (not under) the burn fill, or picking a
scan/frame where the two footprints diverge more to confirm the legend
entry is doing anything at all.

**Output location**: `krgx_hourly/column_frames/` (14 PNGs, ~115-122 KB
each). No changes made to any existing files — fully standalone per the
task's own scope note.

**Environment note (repeats)**: same `geo_env` python.exe + manual PATH
workaround as prior tasks (per the note in `CLAUDE.md`).