"""
2D KDE landscape plot for a Cobaya chain — Tom's "Method 1" diagnostic.

Produces a 3-panel figure on (P1, P2):
  - main:  2D KDE contour heatmap with median × overlay and faint sample dots
  - top:   1D KDE marginal of P1 with median tick
  - right: 1D KDE marginal of P2 with median tick

Run on 3c (Mac mini) data:
    python3 kde_plot.py /tmp/3c_chain afterglow_3c beta_aft c_D

Run on 3b (Ubuntu) data once the chain is pullable:
    python3 kde_plot.py /tmp/3b_chain afterglow_3b beta_aft c_D
"""
import sys, os, glob, json
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

if len(sys.argv) < 5:
    sys.exit("Usage: kde_plot.py <chain_dir> <chain_prefix> <param1> <param2> [--burn 0.3] [--out <png>]")

chain_dir, prefix, p1, p2 = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
burn_frac = 0.3
out_png = os.path.join(chain_dir, f"kde_{prefix}_{p1}_{p2}.png")
for i, a in enumerate(sys.argv[5:]):
    if a == "--burn":  burn_frac = float(sys.argv[5 + i + 1])
    if a == "--out":   out_png   = sys.argv[5 + i + 1]

# --- Load ---
files = sorted(glob.glob(os.path.join(chain_dir, f"{prefix}.[1-9].txt")))
if not files: sys.exit(f"no walker files at {chain_dir}/{prefix}.*.txt")
header = open(files[0]).readline().lstrip("#").split()
if p1 not in header: sys.exit(f"{p1} not in header (have: {header[:5]} ...)")
if p2 not in header: sys.exit(f"{p2} not in header")
iw, i1, i2 = header.index("weight"), header.index(p1), header.index(p2)

raw = [np.loadtxt(f) for f in files]
print(f"# walkers: {len(raw)}, rows: {[r.shape[0] for r in raw]}")
chunks = [r[int(burn_frac * len(r)):] for r in raw]
data = np.vstack(chunks)
w = data[:, iw].astype(int)
print(f"# after {int(burn_frac*100)}% burn-in: {data.shape[0]} unique rows, weight-sum {w.sum()}")

x = np.repeat(data[:, i1], w)
y = np.repeat(data[:, i2], w)

# Subsample for KDE speed and dot overlay
rng = np.random.default_rng(11)
if len(x) > 80_000:
    idx = rng.choice(len(x), 80_000, replace=False); x = x[idx]; y = y[idx]
print(f"# samples after subsample: {len(x)}")

# Robust medians (full sample, not subsampled)
x_full = np.repeat(data[:, i1], w); y_full = np.repeat(data[:, i2], w)
med_x, med_y = float(np.median(x_full)), float(np.median(y_full))
print(f"# medians: {p1}={med_x:.4f}  {p2}={med_y:.4f}")

# --- 2D KDE on a grid ---
xy = np.vstack([x, y])
kde = gaussian_kde(xy, bw_method='scott')
pad_x = 0.05 * (x.max() - x.min()); pad_y = 0.05 * (y.max() - y.min())
gx = np.linspace(x.min() - pad_x, x.max() + pad_x, 160)
gy = np.linspace(y.min() - pad_y, y.max() + pad_y, 160)
GX, GY = np.meshgrid(gx, gy)
GZ = kde(np.vstack([GX.ravel(), GY.ravel()])).reshape(GX.shape)

# Density at median (for in-figure annotation)
dens_at_med = float(kde(np.array([[med_x], [med_y]]))[0])
dens_peak = float(GZ.max())
peak_idx = np.unravel_index(np.argmax(GZ), GZ.shape)
peak_x, peak_y = float(GX[peak_idx]), float(GY[peak_idx])
ratio = dens_at_med / dens_peak

# --- 1D KDEs ---
kx = gaussian_kde(x_full[:80_000] if len(x_full) > 80_000 else x_full)
ky = gaussian_kde(y_full[:80_000] if len(y_full) > 80_000 else y_full)

# --- Figure ---
fig = plt.figure(figsize=(11, 9))
gs = GridSpec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4],
              hspace=0.04, wspace=0.04, left=0.08, right=0.97, top=0.93, bottom=0.08)

ax_main = fig.add_subplot(gs[1, 0])
ax_top  = fig.add_subplot(gs[0, 0], sharex=ax_main)
ax_rt   = fig.add_subplot(gs[1, 1], sharey=ax_main)

# Main: filled contour KDE
levels = np.linspace(0, GZ.max(), 14)[1:]
cf = ax_main.contourf(GX, GY, GZ, levels=levels, cmap='viridis')
ax_main.contour(GX, GY, GZ, levels=levels, colors='white', linewidths=0.3, alpha=0.5)
# Faint sample dots
ax_main.scatter(x[::20], y[::20], s=0.5, c='white', alpha=0.25, edgecolors='none')
# Median ×
ax_main.scatter([med_x], [med_y], marker='X', s=320, c='red', edgecolor='black',
                linewidth=2, zorder=10, label=f'median ({med_x:.3f}, {med_y:.3f})')
# Peak ·
ax_main.scatter([peak_x], [peak_y], marker='*', s=320, c='gold', edgecolor='black',
                linewidth=1.5, zorder=11, label=f'KDE peak ({peak_x:.3f}, {peak_y:.3f})')
ax_main.set_xlabel(p1); ax_main.set_ylabel(p2)
ax_main.legend(loc='upper right', fontsize=10, framealpha=0.9)

# Top marginal: 1D KDE of p1
gx1 = np.linspace(x.min() - pad_x, x.max() + pad_x, 400)
ax_top.fill_between(gx1, kx(gx1), color='steelblue', alpha=0.6)
ax_top.axvline(med_x, color='red', lw=2, ls='--')
ax_top.text(med_x, 0.95 * kx(gx1).max(), f' median {med_x:.3f}', color='red',
            fontsize=9, va='top')
ax_top.set_ylabel("density"); ax_top.tick_params(labelbottom=False)
ax_top.set_title(f"2D KDE landscape — {prefix} — ({p1}, {p2}) — {len(x_full)} weighted samples",
                 fontsize=12, pad=8)

# Right marginal: 1D KDE of p2
gy1 = np.linspace(y.min() - pad_y, y.max() + pad_y, 400)
ax_rt.fill_betweenx(gy1, 0, ky(gy1), color='steelblue', alpha=0.6)
ax_rt.axhline(med_y, color='red', lw=2, ls='--')
ax_rt.text(0.95 * ky(gy1).max(), med_y, f'median\n{med_y:.3f}', color='red',
           fontsize=9, ha='right', va='center')
ax_rt.set_xlabel("density"); ax_rt.tick_params(labelleft=False)

# Diagnostic annotation box
verdict = ("median sits ON the peak" if ratio > 0.9 else
           "median NEAR the peak" if ratio > 0.5 else
           "median in low-density region — possible bimodality" if ratio > 0.2 else
           "median IN THE VALLEY (low density) — bimodality likely")
ax_main.text(0.02, 0.02,
             f"ρ(median) / ρ(peak) = {ratio:.3f}\n→ {verdict}",
             transform=ax_main.transAxes, fontsize=10,
             bbox=dict(facecolor='white', alpha=0.85, edgecolor='gray'),
             va='bottom', ha='left')

# Colorbar
cax = fig.add_axes([0.62, 0.95, 0.30, 0.02])
fig.colorbar(cf, cax=cax, orientation='horizontal', label='density')

fig.savefig(out_png, dpi=130, bbox_inches='tight')
print(f"Wrote {out_png}")

# Also write a small JSON of the headline numbers
out_json = os.path.splitext(out_png)[0] + ".json"
with open(out_json, "w") as f:
    json.dump({
        "chain_prefix": prefix, "p1": p1, "p2": p2,
        "n_weighted_samples": int(w.sum()),
        "median": [med_x, med_y], "kde_peak": [peak_x, peak_y],
        "density_at_median": dens_at_med, "density_at_peak": dens_peak,
        "median_over_peak_ratio": ratio, "verdict": verdict,
    }, f, indent=2)
print(f"Wrote {out_json}")
