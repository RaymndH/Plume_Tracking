# Bug Fire (2026-08-10) Dual-Pol Smoke and Column Tracking (KRGX)

Earlier Claude Code project run in PM/agent mode: a Claude chat wrote tasks, Claude Code executed and appended results to `TASKS.md`. Tasks 01 to 25 complete.

## Setup
- Radar: KRGX (Reno, NV), NEXRAD Level-II dual-pol from S3
- Libs: pyart, numpy, scipy, matplotlib
- Windows conda env `geo_env` (miniforge3). Claude Code's PowerShell runs with `-NoProfile -NonInteractive`, so call geo_env's `python.exe` directly and prepend its dirs to PATH manually

## Decisions
- Smoke mask: RhoHV-only (`RHV_MAX=0.90`); ZDR dropped as too noisy at this SNR
- Fire source anchors: GOES-18 FDCC power-weighted fire-pixel centroids (radar plume centroids caused a circularity bug)
- "Plume 1" reclassified as a secondary lobe of the Bug Fire, not a spot fire (GOES shows no independent fire pixels there)
- Tracker: Bayesian probabilistic occupancy grid (Tasks 16/19); nearest-centroid and IoU trackers failed from lobe swapping
- 20 km/hr speed limit applies to fire-front propagation, not smoke-centroid advection (applying it to smoke caused 8/13 tracking failures)

## Open problem: column tracking (Tasks 21 to 25)
- The rising injection column near the source is not the most intense or polarimetrically distinct feature; aged advected smoke dominates dBZ and RhoHV
- Nearest-blob connected components (Task 25) still land 4 to 7 km from the GOES anchor
- Candidate methods:
  - GOES-seeded marker-controlled watershed (recommended Task 26): GOES anchor as seed marker, grow column boundary
  - Persistent homology (as in the tobac Camp Fire pipeline)
  - ZDR column detection: elevated ZDR in the updraft core above the source

## Outputs
- `krgx_hourly/column_tracker.json`, `hourly_results.json`, `smoke_footprints.npz`
- `grids_*.npz` (14 per-scan grids)
- PNG frames: `centroid_frames/`, `column_frames/`
