import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from datetime import datetime, timezone
from pathlib import Path

# ── load results ──────────────────────────────────────────────────────────
with open('krgx_hourly/hourly_results.json') as f:
    results = json.load(f)

# ── footprint-overlap (IoU) tracking, speed cap kept only for reporting ───
times   = []
x_cents = []
y_cents = []
z_cents = []
z_tops  = []
smoke_fracs = []
scan_labels = []

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

    # track_ok[i] is True if scan i's link to scan i-1 is valid
    # For a speed between scan i and i+1, both endpoints must be valid links
    interval_ok = track_ok[:-1] & track_ok[1:]

    speed_kmph = np.where(
        interval_ok,
        np.sqrt(dx**2 + dy**2) / (dt_min / 60),
        np.nan
    )
    bearing = np.where(interval_ok, bearing, np.nan)

    print(f"\nScan-to-scan drift velocity:")
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

# ── plots ─────────────────────────────────────────────────────────────────
import matplotlib
matplotlib.use('Agg')

fig, axes = plt.subplots(2, 2, figsize=(12, 9))

# 1. X/Y centroid track (plan view)
ax = axes[0, 0]
sc = ax.scatter(x_cents, y_cents, c=times, cmap='plasma', s=60, zorder=3)
ax.plot(x_cents, y_cents, 'k--', lw=0.8, zorder=2)
for i in range(1, len(x_cents)):
    if not track_ok[i]:
        ax.plot([x_cents[i-1], x_cents[i]],
                [y_cents[i-1], y_cents[i]],
                'r--', lw=1.5, zorder=2)
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
