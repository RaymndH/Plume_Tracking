import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
    axes[0].set_xlabel('RhoHV'); axes[0].set_title('Plume 0 - RhoHV')
    axes[0].legend()
if len(zdr_p0):
    axes[1].hist(zdr_p0, bins=30, color='darkorange', edgecolor='k')
    axes[1].axvline(-0.5, color='red', lw=2)
    axes[1].axvline( 1.5, color='red', lw=2, label='ZDR window')
    axes[1].set_xlabel('ZDR (dB)'); axes[1].set_title('Plume 0 - ZDR')
    axes[1].legend()
plt.tight_layout()
plt.savefig('krgx_cache/task04_pol_histograms.png', dpi=120)
print("\nSaved task04_pol_histograms.png")
