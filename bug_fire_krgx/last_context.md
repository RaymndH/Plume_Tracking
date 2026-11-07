read CLAUDE.md first, then claude_initial_prompt.md You are a project manager that will be linked to a Claude Code running on a laptop with all the actual scripts etc.

Both files are already in my context from the uploaded documents. I've fully read both. Here's my situation assessment as your project manager:

KRGX Plume Analysis — Project Manager Briefing
What I Know

I have full context from both files. The project is well-defined:

Current state: Working single-scan pipeline (krgx_plume.py) + hourly batch (krgx_hourly.py) using histogramdd gridding, watershed hotspot detection, skeleton extraction, and Plotly 3D viz. GOES-18 fire detection and NIFC perimeters are integrated.
Target event: Bug & Stallion Fire, 2026-08-09 ~23:59Z, KRGX (Reno)
Known gap: Pipeline is currently reflectivity + spectrum width only — no polarimetric fields (ρHV, ZDR, ΦDP) yet, which the initial prompt identifies as the highest-priority next step
The Strategic Gap

The current algorithm is essentially a Z + SW morphological tracker. The research roadmap (from claude_initial_prompt.md) calls for a full polarimetric masking pipeline before any ML/optical flow work makes sense. Without ρHV and ZDR, you can't reliably separate smoke from precipitation, ground clutter, or biological echoes.

What I Need from You to Get Started

Before I start directing Claude Code to write or modify anything, I need to resolve a few things:

What does KRGX actually have? NEXRAD WSR-88D dual-pol has been standard since ~2013, so KRGX should have ZDR, ρHV, ΦDP, and Vr. But are you currently loading those fields in krgx_plume.py, or only Z and SW? I need to see the file to know what's already there.
What's the immediate goal for this session? Options:
Add polarimetric fields (ρHV, ZDR) to the existing histogramdd grid
Build the fuzzy-logic smoke mask (Path 1 from the roadmap)
Extend krgx_hourly.py to output a time series of blob stats
Something else you have in mind today
Can you paste or have Claude Code read krgx_plume.py? That's the central file — I need to see the current gate extraction loop to know exactly where to slot in polarimetric fields without breaking the histogramdd pipeline.

Once I have those answers (especially #3), I can give Claude Code precise, file-aware instructions rather than general guidance.

krgx_plume.py
PY

Good — I have the full picture now. Here's a complete audit:

krgx_plume.py — Code Audit & Gap Analysis
What the Code Does Well

The pipeline is clean and consistent. The histogramdd approach is the right call — the CLAUDE.md decision to abandon grid_from_radars was correct. The SW isosurface rendering is solid, the AGL→MSL conversion is handled correctly throughout, and the gate extraction loop is well-structured.

Critical Gap: Polarimetric Fields Are Not Loaded

This is the central issue. Line 161:

python
radar = pyart.io.read_nexrad_archive(
    io.BytesIO(raw), include_fields=['reflectivity', 'spectrum_width']
)

ZDR, ρHV, ΦDP, and Vr are explicitly excluded. KRGX is a dual-pol WSR-88D — all four fields are in the Level-II archive. They're just never requested.

This is the only thing blocking Path 1 (polarimetric smoke mask). Everything else — the gate extraction loop, the histogramdd gridding, the Gaussian smooth — is already structured to accept additional fields. Adding polarimetric fields is a surgical change, not a rebuild.

Secondary Issues
Issue	Location	Severity
sw variable name collision — loop variable sw (sweep index) overwritten by sw (spectrum width field data)	Line 183–186	Bug — sw_field['data'][sl] uses the wrong loop variable name
radar deleted before PCA/DBSCAN (line 204) — good practice, but note it	Line 204	OK
Terrain meshgrid uses (nx, ny) shape but terrain is (ny, nx) — line 323 YY_t, XX_t = np.meshgrid(y_km, x_km) gives (ny, nx) arrays, but terrain is also (ny, nx), and then terrain_km.T is used in a commented-out trace but terrain_km (no transpose) in the live trace (line 463)	Lines 323, 463	Possible shape mismatch — needs verification
binary_3d = dbz_smooth > 15 hardcodes 15 dBZ, ignoring DBZ_THRESHOLD = 5 config	Line 278	Minor inconsistency
The Variable Name Bug (Important)

In the gate extraction loop at lines 176–202:

python
for sw in keep:           # ← sw = sweep index (int)
    sl = radar.get_slice(sw)
    ...
    sw_field = radar.fields.get('spectrum_width')
    if sw_field is not None:
        sw = sw_field['data'][sl]   # ← sw is now overwritten with array!

The sweep index sw gets clobbered by the spectrum width array on the first iteration. This means radar.get_slice(sw) would fail on iteration 2 because sw is now a numpy array. This is either masked by the fact that only one sweep survives the time window (single scan), or it's silently broken for multi-sweep runs.

What I'd Direct Claude Code to Do Next

Here's the precise task list, in order:

Task 1 — Fix the variable name collision (5 min)
Rename loop variable sw → isw throughout the loop.

Task 2 — Add polarimetric fields to include_fields and gate extraction
Extend include_fields to ['reflectivity', 'spectrum_width', 'differential_reflectivity', 'cross_correlation_ratio', 'differential_phase'], add corresponding extraction in the gate loop (same valid mask), and add histogramdd passes for each field.

Task 3 — Build a smoke_mask_3d boolean array
After gridding: ρHV < 0.75 AND ZDR > 0.5 dB AND dBZ > DBZ_THRESHOLD → this is your polarimetric smoke discriminant, grounded in Melnikov 2008/2020.

Task 4 — Add the smoke mask isosurface to the Plotly figure
One new go.Isosurface or go.Volume trace colored by ρHV, with the smoke mask applied.

Do you want me to write up the exact code diff for Claude Code to execute (Tasks 1–3 first)? Or do you want to discuss the polarimetric threshold choices before touching the code?

Good question — this is the right thing to nail down before we start generating code specs.

Best Workflow: PM (me) ↔ Claude Code
The Core Problem

Claude Code on your laptop and me here are separate instances with no shared memory. Every handoff has to be explicit text — either pasted into Claude Code's terminal, or written to a file it can read.

Recommended Setup: TASKS.md as the Shared State File

Create C:\Users\raymondharding\Research\Runs\TASKS.md. This file is the contract between us.

My role: I write precise, scoped task specs into this chat. You copy them into TASKS.md (or I write them directly if you paste the file back to me).

Claude Code's role: Reads TASKS.md, executes the current task, appends results/questions to the bottom, stops.

Your role: Paste Claude Code's output/questions back here. I update the spec or write the next task.

File Structure I'd Suggest
Runs/
├── CLAUDE.md              ← already exists (project context)
├── TASKS.md               ← NEW: current task queue + results log
├── krgx_plume.py          ← existing
├── krgx_hourly.py         ← existing
└── krgx_modules/          ← NEW: one file per module
    ├── pol_classifier.py  ← Module 1
    ├── plume_geometry.py  ← Module 2/3
    ├── goes_fusion.py     ← Module 4
    └── fei_tracker.py     ← Module 5
How a Cycle Works
Me (here)          →  writes TASK spec to chat
You                →  paste into TASKS.md, tell Claude Code "read TASKS.md and execute Task N"
Claude Code        →  writes code, appends output/errors/questions to TASKS.md
You                →  paste TASKS.md back here (or just paste the results section)
Me (here)          →  reviews, writes next task or fixes spec

The key discipline: Claude Code appends, never overwrites TASKS.md. The file accumulates a log of what was done and what was found. That log is what I read to write the next task.

First TASKS.md — Ready to Paste

Here's your starting file. Copy this exactly into TASKS.md:

markdown
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
(Claude Code pastes output here)

That's it. Paste that file, run Task 01, and paste the results section back here. I'll write Task 02 based on what fields are actually available.

Want me to also draft Task 02 speculatively (assuming full dual-pol) so we can move faster once Task 01 confirms?