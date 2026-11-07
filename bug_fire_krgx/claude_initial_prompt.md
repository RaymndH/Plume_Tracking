# claude_initial_prompt.md
## Wildfire Radar Plume Detection — Project Handoff Document
*Generated from initial research conversation, 2026-09-18*

---

## Project Context

This document captures the full context of the initial research and planning session for a wildfire smoke plume detection system using weather radar data. It is intended as a handoff to an agent managing this project going forward.

The researcher (Raymond Harding) has existing Python code and Jupyter notebooks in `C:\Users\raymondharding\Research\Runs`. Those files were not uploaded in this session — they should be read directly by Claude Code or pasted in by the user.

The core goal is to **move away from explicit algorithmic Python toward a data science / machine vision / statistical pipeline** for identifying wildfire smoke plumes in weather radar data, characterizing their structure, and estimating fire behavior.

---

## What Was Discussed

### 1. The Research Question

Given weather radar data (NEXRAD WSR-88D or similar), how do you:
- Correctly identify which radar signals belong to smoke plumes (vs. precipitation, ground clutter, biological echoes)
- Characterize plume structure (area, orientation, vertical extent)
- Estimate plume rate of travel and direction
- Identify possible spotting structure (secondary downwind echoes)
- Filter noise using Digital Elevation Models (DEM) for beam blockage
- Incorporate fuel maps and slope to estimate Rate of Spread (ROS)
- Use wind data to infer plume direction

### 2. Key Analogies: Neuroscience / fMRI

Radar volume scans over time are structurally similar to fMRI BOLD data: a 3D spatial volume sampled at discrete time steps, with signals of interest embedded in noise. Methods from neuroimaging transfer directly:

| Neuroscience / fMRI | Radar Wildfire |
|---|---|
| BOLD signal intensity | Radar reflectivity Z |
| Activated voxels vs baseline | Smoke echoes vs clear air / precipitation |
| Optical flow of traveling waves | Plume advection velocity |
| Source/sink identification | Fire origin vs downwind dispersion |
| Spatial smoothing / GLM masking | Polarimetric thresholding + DEM blockage |
| Temporal SNR (BOLD response curve) | FEI temporal evolution (rise / plateau / decay) |
| Multi-voxel pattern analysis | Multi-variable polarimetric feature space |
| ConvRNN for brain dynamics | ConvLSTM for plume forecasting |

Key neuroscience methods that apply:
- **Optical flow** (Horn-Schunck, CLG, Farneback): used in mesoscale brain imaging to detect traveling waves and source/sink locations. Apply to successive radar PPI frames to extract plume motion vectors.
- **Spatiotemporal pattern detection**: pre-filter (Gaussian/median spatial smoothing), normalize across elevation angles, then apply optical flow. The plume is the "traveling wave."
- **ConvLSTM / spatiotemporal deep learning**: multi-layer attention modules address long-term feature dependencies in sequential echo images (ConvRNNs outperform traditional optical flow but need attention for intensity discrimination).

References:
- Optical flow in neuroscience: Afrashteh et al. 2017 (NeuroImage); Townsend & Gong 2018 (PLoS Comp Bio)
- ConvRNN for radar: MLAM paper (HKO-7 / HMB datasets)

---

## Polarimetric Radar Signatures of Smoke Plumes

This is the most important anchor for plume masking. Research consensus (Melnikov et al. 2008, 2009, 2020; Brown et al. 2021; Wu et al. 2019):

### Copolar Correlation Coefficient (ρHV / rho_hv) — PRIMARY MASK
- S-band observations: ρHV in wildfire plumes mostly **< 0.4** (Melnikov 2008, Oklahoma WSR-88D)
- X-band (southern China): ρHV **0.55–0.75**, average 0.61 (higher at X-band due to Mie scattering)
- Caused by highly nonspherical, randomly-oriented ash particles with wide size distributions
- **Use:** Threshold ρHV < 0.6–0.75 (band-dependent) to isolate non-meteorological echoes likely to be smoke/ash
- Note: ρHV in smoke is statistically biased due to low correlations and number of samples M — account for SNR-dependent bias

### Differential Reflectivity (ZDR) — SECONDARY DISCRIMINANT
- Smoke plumes: ZDR **> 0**, typically ~2 dB (S-band); 2.28–5.08 dB average 3.99 dB (X-band China)
- Indicates oblate/horizontally-oriented scatterers (ash particles, debris)
- Positive ZDR + low ρHV = strong discriminant vs. many other non-met echoes
- **Use:** ZDR > 0.5 dB as secondary gate alongside ρHV threshold

### Spectrum Width (σv) — TURBULENCE INDICATOR
- Consistent spectrum widths in plumes ~2 m/s (S-band)
- Captures velocity dispersion within the beam — turbulent pyrogenic convection creates wide σv
- **Use:** High σv collocated with low ρHV + positive ZDR → active/turbulent plume core and probable spotting zones

### Differential Phase (ΦDP)
- Southern China X-band: ΦDP **40.27°–69.27°**, average 56.66°
- Large ΦDP values help separate smoke from precipitation
- **Use:** Additional soft discriminant in multi-variable fuzzy classifier

### Temporal Elevation Lag — PLUME GROWTH DIAGNOSTIC
- Changes in echo intensity at higher elevation scans **lag behind** lower elevation scans
- Change in wildfire area lags behind variation in echo intensity
- Lower elevation detects plume first (~5 min after ignition at 1.8°), then higher tilts fill in as plume rises
- **Use:** This lag structure is a diagnostic signature of active plume growth vs. advected/decaying smoke

### Fire Echo Intensity (FEI) Index
- Composite index combining ZH, ZDR, and ρHV
- Temporal evolution: rapid rise during early plume ascent (particle concentration increases) → gradual decline as buoyancy weakens and particles settle
- **Use:** Engineered feature for ML classifier; temporal shape is a fire activity proxy

---

## DEM-Based Beam Blockage

Mountainous terrain causes partial or complete radar beam blockage — critical for the Sierra Nevada / western US context.

### Method
1. For each radar site, compute a **2D blockage fraction map** (azimuth × range) at each tilt angle using:
   - SRTM or 3DEP DEM data
   - 4/3 Earth radius beam propagation model
   - For each range bin: compute upper/lower beam edge elevation, compare to terrain elevation
   - If terrain exceeds lower beam edge → beam is blocked
2. Mask out range-azimuth cells where blockage fraction exceeds threshold (e.g., >50% power loss)
3. Stack across tilt angles → **3D "unblocked" volume** — only use voxels with reliable coverage
4. This blockage map is **static** (precompute once per radar site)

### References
- GIS-based DEM methodology: Krajewski et al. 2005 (Computers & Geosciences) — NWS Beam Blockage Algorithm (BBA) basis
- 4/3 Earth radius model standard in NWS NEXRAD Product Generation System

### Fuel Canopy as Additional Blockage
- Dense forest canopy attenuates signal before it reaches surface
- Fuel map canopy cover × biomass → soft prior for additional attenuation
- Apply as a weight layer on top of terrain blockage (not a hard mask)

---

## Fuel Maps, Rate of Spread, and Slope

### Rothermel ROS Model
The Rothermel fire spread model (1972) defines Rate of Spread as:

```
R = R₀ × (1 + Φw + Φs)
```

Where:
- R₀ = base spread rate (no wind, no slope)
- Φw = wind correction factor (function of mid-flame wind speed)
- Φs = slope correction factor

Wind and slope are the two dominant variables. Their effects are combined by vector addition.

### Using ROS as a Spatial Prior
- DEM slope + wind direction → expected plume elongation axis (major axis aligns with wind-slope combined vector)
- High-ROS zones (steep slope + aligned wind + dry fuel) → expected spotting and turbulence → cross-validate with high σv regions
- Fuel load (LANDFIRE) → estimate scatterer density (more fuel = more pyrometeors = higher Z for same fire intensity)
- Use ROS map to weight plume probability higher in physically plausible zones

### Data Sources
- **DEM**: SRTM 30m or USGS 3DEP 10m/1m
- **Fuel models**: LANDFIRE (Scott & Burgan 40 fuel models, canopy cover, canopy height, canopy bulk density)
- **Wind**: HRRR, RAP, or surface observation networks
- **Slope**: computed from DEM (aspect + gradient magnitude)

---

## Machine Learning / Machine Vision Pipeline

### Fuzzy Logic Classifier (Interpretable Baseline)
Following Wu et al. 2019 (southern China X-band paper):
- Input features: Z, ZDR, ρHV, σv, ΦDP, temporal Z gradient
- Membership functions tuned to published smoke plume distributions (Gaussian fits to Melnikov 2008/2020 histograms)
- Three-class output: **smoke plume / precipitation / ground clutter**
- Weight matrix optimized per radar band (S vs. X vs. C)
- Fast, interpretable, no labeled training data needed initially

### Feature Engineering Per Voxel
For each range-azimuth-elevation voxel:
- Z, ZDR, ρHV, ΦDP, σv, Vr (radial velocity)
- Spatial gradient of Z and ρHV (plume edge detection)
- Temporal derivative of Z (FEI growth rate)
- DEM blockage fraction (mask/weight)
- Wind alignment score: dot product of local Vr with NWP wind vector at that altitude
- FEI composite index

### Blob Detection / Connected Components
After polarimetric masking:
1. Binary smoke mask from threshold gates
2. Connected component labeling (scipy.ndimage or skimage)
3. Per-blob features: area, centroid, bounding ellipse, major/minor axis (orientation = plume direction), eccentricity, aspect ratio
4. Track blobs across time: Hungarian algorithm or Kalman filter for assignment
5. FEI time series per blob

### Optical Flow for Velocity
- Lucas-Kanade or Farneback dense optical flow on sequential masked Z fields
- Per-pixel velocity vectors → average over blob → plume advection speed + direction
- Divergence/convergence fields → identify active growth zones
- Compare blob velocity to surface wind → estimate pyrogenic uplift (deviation from surface flow = strong pyroconvection)

### ConvLSTM for Plume Forecasting
- Input: sequence of 3D polarimetric volumes (Z, ZDR, ρHV, σv) — typically 5–10 frames
- Output: predicted plume mask at T+N minutes
- Train on NEXRAD Level 2 archive data for known fire events (historical California fires have extensive coverage)
- Attention mechanism to weight turbulent core

### Spotting Structure Detection
Secondary blob detection pipeline:
1. Identify isolated blobs downwind from main plume
2. Feature gates: aligned with wind vector, separated by > N km from main plume, FEI growing, σv elevated
3. Cross-validate with fuel availability (LANDFIRE — spotting requires receptive fuel)
4. Output: spotting probability map

---

## Key References

### Chinese X-Band Radar Paper (PRIMARY)
**Wu et al. 2019 — "Monitoring wildfire using high-resolution compact X-band dual-polarization radar: A case study in southern China"**
- *Remote Sensing of Environment / ScienceDirect*
- Developed fuzzy-logic Wildfire Identification Algorithm (WIA) using Z, ZDR, ρHV, std(Vr)
- Key findings: ρHV 0.55–0.75, ZDR 2.28–5.08 dB, ΦDP 40.27–69.27°
- Temporal-elevation lag diagnostic for plume growth
- Demonstrated X-band compact radar as viable alternative to long-range S-band

### Melnikov et al. 2008
- "Radar polarimetric signatures of fire plumes in Oklahoma" — *Geophysical Research Letters*
- S-band WSR-88D observations: ρHV < 0.4, positive ZDR
- Established the polarimetric discriminant values that are now standard

### Melnikov et al. 2009
- "Polarimetric radar properties of smoke plumes: A model" — *JGR Atmospheres*
- Physical model of scatterer shapes (needle-like at 10 cm wavelength)
- Statistical bias in ρHV estimation for smoke (important for threshold setting)

### Melnikov et al. 2020 ("Of Fire and Smoke Plumes, Polarimetric Radar Characteristics")
- *Atmosphere* (MDPI)
- Comprehensive survey across grass, brush, forest fires
- Notable: ZDR ~2 dB, ρHV ~0.6, ΦDP large, σv ~2 m/s, Vr ~20 m/s smooth field
- Comparison S-band vs. C-band (TDWR)

### Brown et al. 2021 (Ka-Band Mobile Radar)
- "Mobile Ka-Band Polarimetric Doppler Radar Observations of Wildfire Smoke Plumes"
- *Monthly Weather Review*
- Fine-scale plume kinematics, Z -15 to 20 dBZ, Vr 0–16 m/s
- ρHV < 0.8 in plume; oblate-shaped scatterers

### Du et al. 2025
- "A Technology of Forest Fire Smoke Detection Using Dual-Polarization Weather Radar"
- *Forests* (MDPI)
- FEI index temporal evolution; recommends adding ΦDP, Vr, σv to future work

### DEM Blockage
- Krajewski et al. 2005 — GIS-based DEM beam blockage methodology (Computers & Geosciences)
- NWS Beam Blockage Algorithm (BBA) — NEXRAD standard, 4/3 Earth radius model

### Rothermel / Fire Spread
- Rothermel 1972 (original ROS model)
- Andrews 2012 — "Modeling Wind Adjustment Factor and Midflame Wind Speed"
- WRF-Fire paper (Mandel et al.) — coupled atmosphere-wildland fire model implementing Rothermel

---

## Paths Forward (Prioritized)

### Path 1 — Polarimetric Masking + Statistical Blob Characterization (Weeks 1–4)
Quickest win. Replaces most of the current pixel-by-pixel Python with a statistically grounded filter.

1. Load NEXRAD Level 2 via `pyart` or `wradlib`
2. Apply DEM blockage mask (precomputed, static — use `pyproj` + SRTM)
3. Threshold: ρHV < 0.7, ZDR > 0.5 dB, Z > 5 dBZ → binary smoke mask
4. Connected component labeling (`scipy.ndimage.label` or `skimage.measure.label`)
5. Extract per-blob: area, centroid, orientation (major axis = plume direction), eccentricity, FEI
6. Output: structured DataFrame of blob properties per scan time

### Path 2 — Optical Flow Velocity Estimation (Weeks 3–6)
Apply to successive masked Z fields (OpenCV `calcOpticalFlowFarneback` or `scikit-image`):
- Per-frame plume velocity vectors
- Plume advection speed and direction
- Divergence/convergence → active growth zones
- Deviation from NWP wind → pyrogenic uplift estimate

### Path 3 — DEM + Fuel Integration (Weeks 4–8)
1. Compute 2D blockage fraction maps per tilt: `pyproj` beam propagation + `rasterio` DEM read
2. Ingest LANDFIRE fuel model + canopy cover → soft attenuation layer
3. Compute per-pixel Rothermel ROS proxy (DEM slope + NWP wind + fuel model)
4. Use ROS map as spatial prior weighting plume probability
5. Cross-validate plume orientation against ROS gradient direction

### Path 4 — Fuzzy Logic Classifier (Weeks 6–10)
Implement WIA-style fuzzy classifier (Wu et al. 2019):
- Features: Z, ZDR, ρHV, σv, ΦDP, temporal ΔZ
- Membership functions from published Gaussian fits to plume histograms
- Three-class output: smoke / precip / clutter
- Interpretable, fast, no labeled training data needed initially

### Path 5 — ConvLSTM Plume Tracker (Weeks 10–20)
Once labeled dataset is built from Paths 1–3:
- Train ConvLSTM on NEXRAD Level 2 archive (known California fire events)
- Input: 5–10 frame sequences of (Z, ZDR, ρHV, σv) at multiple tilts → CAPPI or constant-altitude grids
- Output: plume mask + FEI growth prediction
- Attention mechanism for turbulent core / spotting zone weighting

### Path 6 — Spotting Structure Detection
Secondary-blob pipeline (can run alongside Paths 1–2):
- Candidate spots: isolated blobs downwind from main plume
- Feature gates: wind alignment, distance from main plume, FEI growing, σv elevated
- Cross-validate with LANDFIRE fuel availability
- Output: spotting probability map

---

## Open Questions (Need Files to Answer)

1. **What does the current Python algorithm do?** — need `.py` files from `C:\Users\raymondharding\Research\Runs`
2. **Which radar network?** NEXRAD WSR-88D? State network? X-band? — determines available polarimetric variables and resolution
3. **What geographic area?** — determines relevant radar sites, DEM data source, fuel model coverage
4. **What is the current output format?** — what does the algorithm produce today?
5. **What are the known failure modes?** — when does the current algorithm miss plumes or produce false positives?
6. **Are there labeled examples?** — confirmed plume events with timestamps and locations for validation/training?

---

## Agent Instructions

When the user provides the files from `C:\Users\raymondharding\Research\Runs`:

1. Read `claude.md` first for any additional context Raymond has added
2. Read `.ipynb` notebooks to understand the current workflow, data loading, and visualization
3. Read `.py` files to understand the current algorithm logic
4. Map the current algorithm to the paths above — identify which parts already exist and where the gaps are
5. Prioritize Path 1 (polarimetric masking + blob characterization) as the first concrete deliverable
6. All new code should be data/ML-centered, not explicit rule-based Python where avoidable
7. Use `pyart` for NEXRAD Level 2 I/O, `rasterio` for DEM, `scipy`/`skimage` for blob detection, `OpenCV` for optical flow
8. Keep the DEM blockage map as a precomputed static asset (recompute only when radar site changes)
9. The FEI index (Z + ZDR + ρHV composite) should be a first-class output at every scan time
10. Document all threshold choices with citations to the polarimetric radar literature above
