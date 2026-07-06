"""
Bimodality + median diagnostics for a Cobaya chain directory.

Designed to be droppable on the Ubuntu workstation (3a / 3b) or the Mac mini
(3c) and produce a single JSON the user can paste back for assembly.

Usage:
    python3 analyze_bimodality.py <chain_dir> <chain_prefix> [<param_a> <param_b>] [--burn 0.3]

Examples (paths illustrative):
    # 3b (Ubuntu) — the bimodality Tom is asking about
    python3 analyze_bimodality.py ~/runs/phase3b_noncrossing afterglow_3b beta c_D

    # 3a (Ubuntu) — baseline sanity-check, no beta
    python3 analyze_bimodality.py ~/runs/phase3a_baseline   afterglow_3a c_D w_X

    # 3c (Mac mini) — early-stage sample, B1-narrow prior
    python3 analyze_bimodality.py ~/runs/phase3c_b1narrow   afterglow_3c beta c_D

Outputs:
    diagnostics_<chain_prefix>.json   — per-target univariate + joint stats
    bimodality_<chain_prefix>.png     — KDE marginals + 2D scatter

Decision rule used (call something genuinely bimodal):
    1. ΔBIC(k=1) − ΔBIC(k=2) > 10   (strong preference for mixture)
    2. Dip-test p < 0.01            (rejects unimodality directly)
    3. Density valley exists between the two GMM means on the line connecting
       them, with valley/peak ratio < 0.7  (rules out skewed unimodal mis-fit)
    4. Per-walker marginal medians DISAGREE by more than 2× the within-walker
       std/sqrt(N_per_walker)       (rules out 1-walker outlier squeezing into
                                     a tail that looks like a second mode)

All four must hold to call it "real bimodal." Three of four → "ambiguous,
keep watching." Fewer → "unimodal."
"""
import sys, os, glob, json, argparse
import numpy as np
from scipy.stats import gaussian_kde
from sklearn.mixture import GaussianMixture
import diptest
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("chain_dir")
ap.add_argument("chain_prefix")
ap.add_argument("params", nargs="*", help="Parameter names to analyze in 2D (default: try beta, c_D, then w_X)")
ap.add_argument("--burn", type=float, default=0.3, help="Fraction of each walker dropped as burn-in")
ap.add_argument("--max-samples-2d", type=int, default=120_000, help="Subsample joint analysis above this")
ap.add_argument("--out-dir", default=None)
args = ap.parse_args()

out_dir = args.out_dir or args.chain_dir
os.makedirs(out_dir, exist_ok=True)

# --- 1. Load ---
pattern = os.path.join(args.chain_dir, f"{args.chain_prefix}.[1-9].txt")
files = sorted(glob.glob(pattern))
if not files:
    sys.exit(f"No chain files match {pattern}")
print(f"Found {len(files)} walker files")

header = open(files[0]).readline().lstrip("#").split()
print(f"Columns ({len(header)}): first 12 = {header[:12]} ...")

raw = [np.loadtxt(f) for f in files]
print(f"Rows per walker: {[r.shape[0] for r in raw]}")

# Burn-in
chunks = [r[int(args.burn * len(r)):] for r in raw]
walker_id = np.concatenate([np.full(len(c), i) for i, c in enumerate(chunks)])
data = np.vstack(chunks)
w = data[:, header.index("weight")].astype(int)
print(f"After burn-in: {data.shape[0]} unique rows, weight-sum {w.sum()}")

# Pick joint params
if args.params:
    p1, p2 = args.params[:2]
else:
    candidates = [("beta", "c_D"), ("c_D", "w_X"), ("log10_c_D", "w_X")]
    for a, b in candidates:
        if a in header and b in header:
            p1, p2 = a, b; break
print(f"Joint axes: ({p1}, {p2})")

# Targets for 1D battery
targets = [p for p in [p1, p2, "w_X", "sigma8", "H0", "Omega_m", "c_D"] if p in header]
targets = list(dict.fromkeys(targets))  # de-dup, preserve order

def column(name): return data[:, header.index(name)]
def expand(name, n_max=200_000):
    x = np.repeat(column(name), w)
    if n_max and len(x) > n_max:
        rng = np.random.default_rng(42); x = rng.choice(x, n_max, replace=False)
    return x
def expand_walker_id():
    return np.repeat(walker_id, w)

out = {
    "chain_dir": args.chain_dir, "chain_prefix": args.chain_prefix,
    "n_walkers": len(files), "n_weighted_samples": int(w.sum()),
    "burn_frac": args.burn, "joint_axes": [p1, p2],
    "univariate": {}, "joint": {}, "per_walker_medians": {},
}

# --- 2. Univariate battery ---
print()
print(f"{'param':10s} {'median':>10s} {'mode':>10s} {'ρmed/ρmode':>11s} "
      f"{'dip-p':>8s} {'ΔBIC':>10s}  verdict")
for name in targets:
    x = expand(name)
    if len(x) < 50:
        continue
    kde = gaussian_kde(x, bw_method='scott')
    median = float(np.median(x))
    grid = np.linspace(x.min(), x.max(), 1500)
    dens = kde(grid)
    mode = float(grid[dens.argmax()])
    dens_at_median = float(kde(median)[0])
    dens_at_mode = float(dens.max())
    ratio = dens_at_median / dens_at_mode
    dip_x = x[:60_000] if len(x) > 60_000 else x
    dip_stat, dip_p = diptest.diptest(dip_x)
    Xr = x.reshape(-1, 1)
    g1 = GaussianMixture(1, random_state=0).fit(Xr)
    g2 = GaussianMixture(2, random_state=0).fit(Xr)
    bic = g1.bic(Xr) - g2.bic(Xr)

    # Density-valley test (between the two GMM means, 1D)
    mu = np.sort(g2.means_.flatten())
    line = np.linspace(mu[0], mu[1], 100).reshape(-1, 1)
    dens_line = np.exp(g2.score_samples(line))
    if mu[1] > mu[0]:
        valley_min = float(dens_line.min())
        peaks_min = float(min(dens_line[0], dens_line[-1]))
        valley_ratio = valley_min / peaks_min if peaks_min > 0 else 1.0
    else:
        valley_ratio = 1.0

    real_bimodal = (bic > 10 and dip_p < 0.01 and valley_ratio < 0.7)
    verdict = "BIMODAL" if real_bimodal else (
        "skewed/ambiguous" if bic > 10 else "unimodal"
    )
    out["univariate"][name] = dict(
        median=median, mode=mode, dens_at_median=dens_at_median,
        dens_at_mode=dens_at_mode, density_ratio_median_over_mode=ratio,
        dip_p=float(dip_p), delta_bic_k1_minus_k2=float(bic),
        valley_ratio=float(valley_ratio),
        gmm_k2_means=mu.tolist(), gmm_k2_weights=g2.weights_.tolist(),
        verdict=verdict,
    )
    print(f"{name:10s} {median:10.4f} {mode:10.4f} {ratio:11.2f} "
          f"{dip_p:8.3g} {bic:10.1f}  {verdict}")

# --- 3. Per-walker medians ---
print()
print("Per-walker marginal medians (disagreement => walker stuck in a mode):")
wid_x = expand_walker_id()
for name in targets:
    x = expand(name)
    per_w = []
    for k in range(len(files)):
        mask = wid_x == k
        if mask.sum() > 0:
            per_w.append(float(np.median(x[mask])))
    spread = max(per_w) - min(per_w)
    out["per_walker_medians"][name] = {"medians": per_w, "spread": float(spread)}
    print(f"   {name:10s} per-walker = {[f'{v:.4f}' for v in per_w]}  spread={spread:.4f}")

# --- 4. Joint 2D analysis on (p1, p2) ---
print()
xa, xb = expand(p1, n_max=None), expand(p2, n_max=None)
xy = np.column_stack([xa, xb])
if len(xy) > args.max_samples_2d:
    rng = np.random.default_rng(7)
    idx = rng.choice(len(xy), args.max_samples_2d, replace=False)
    xy = xy[idx]

g1 = GaussianMixture(1, random_state=0, covariance_type='full').fit(xy)
g2 = GaussianMixture(2, random_state=0, covariance_type='full').fit(xy)
joint_bic = g1.bic(xy) - g2.bic(xy)
joint_median = np.median(xy, axis=0)
mu1, mu2 = g2.means_
d_at_mu1 = float(np.exp(g2.score_samples(mu1.reshape(1,-1)))[0])
d_at_mu2 = float(np.exp(g2.score_samples(mu2.reshape(1,-1)))[0])
d_at_med = float(np.exp(g2.score_samples(joint_median.reshape(1,-1)))[0])
# 2D valley test along line connecting the two GMM means
line2d = np.linspace(mu1, mu2, 100)
dens_line = np.exp(g2.score_samples(line2d))
valley_min = float(dens_line.min())
peaks_min = float(min(dens_line[0], dens_line[-1]))
valley_ratio_2d = valley_min / peaks_min if peaks_min > 0 else 1.0

# Per-walker mode-occupancy fractions (the strongest bimodality test)
# Build full-weighted joint sample with walker IDs
xa_full = np.repeat(column(p1), w)
xb_full = np.repeat(column(p2), w)
wid_full = np.repeat(walker_id, w)
xy_full = np.column_stack([xa_full, xb_full])
labels = g2.predict(xy_full)
occupancy = {}
for k in range(len(files)):
    mask = wid_full == k
    if mask.sum() == 0: continue
    frac_B = float((labels[mask] == 1).mean())  # cluster index 1 = mode B
    occupancy[f"walker_{k}"] = frac_B
occ_vals = list(occupancy.values())
occ_spread = max(occ_vals) - min(occ_vals) if occ_vals else 0.0
# If all walkers have similar fraction in B (spread < 0.10), the second cluster
# is just a tail being modeled, not a real basin with stuck walkers.
mixed_well_across_walkers = occ_spread < 0.10
# But the joint bimodality may still be REAL even with good mixing — modes can
# both be visited by every walker. So we keep the geometric test (valley_ratio)
# as the primary criterion, and report walker-mixing as an additional flag.
real_bimodal_joint = (joint_bic > 10 and valley_ratio_2d < 0.7
                      and d_at_med / max(d_at_mu1, d_at_mu2) < 0.5)
likely_tail_not_basin = (real_bimodal_joint and mixed_well_across_walkers
                         and min(g2.weights_) < 0.15)

out["joint"] = {
    "delta_bic_k1_minus_k2": float(joint_bic),
    "k2_means": [mu1.tolist(), mu2.tolist()],
    "k2_weights": g2.weights_.tolist(),
    "k2_covariances": g2.covariances_.tolist(),
    "density_at_meanA": d_at_mu1, "density_at_meanB": d_at_mu2,
    "density_at_pooled_median": d_at_med,
    "median_over_peak_ratio": d_at_med / max(d_at_mu1, d_at_mu2),
    "valley_ratio_2d": valley_ratio_2d,
    "real_bimodal_joint": bool(real_bimodal_joint),
    "per_walker_fraction_in_modeB": occupancy,
    "occupancy_spread": float(occ_spread),
    "well_mixed_across_walkers": bool(mixed_well_across_walkers),
    "likely_tail_not_basin": bool(likely_tail_not_basin),
    "pooled_median": joint_median.tolist(),
}
print(f"JOINT ({p1}, {p2}):  ΔBIC={joint_bic:.1f}  "
      f"valley/peak={valley_ratio_2d:.2f}  "
      f"ρ(med)/max(ρ(mean))={d_at_med/max(d_at_mu1,d_at_mu2):.2f}  "
      f"=>  {'TRULY BIMODAL' if real_bimodal_joint else 'NOT bimodal'}")
print(f"   mean A: ({mu1[0]:.3f}, {mu1[1]:.3f})   weight {g2.weights_[0]:.2f}")
print(f"   mean B: ({mu2[0]:.3f}, {mu2[1]:.3f})   weight {g2.weights_[1]:.2f}")
print(f"   pooled median: ({joint_median[0]:.3f}, {joint_median[1]:.3f})  ← Tom's headline")
print(f"   per-walker fraction in mode B: {[f'{v:.2f}' for v in occ_vals]} "
      f"(spread {occ_spread:.2f})  "
      f"=> {'well-mixed across walkers' if mixed_well_across_walkers else 'WALKER DISAGREEMENT (real bimodality with stuck walker)'}")
if likely_tail_not_basin:
    print(f"   NOTE: small mode weight ({min(g2.weights_):.2f}) + good walker mixing "
          f"→ mode B may be a tail being modeled, not a true basin")

# --- 5. Save & plot ---
out_json = os.path.join(out_dir, f"diagnostics_{args.chain_prefix}.json")
with open(out_json, "w") as f:
    json.dump(out, f, indent=2, default=float)
print(f"\nWrote {out_json}")

# 2x2 plot: KDE(p1), KDE(p2), joint scatter, line-density profile
fig, ax = plt.subplots(2, 2, figsize=(11, 9))
for axi, name in [(ax[0,0], p1), (ax[0,1], p2)]:
    x = expand(name)
    kde = gaussian_kde(x)
    grid = np.linspace(x.min(), x.max(), 600)
    axi.plot(grid, kde(grid), 'k-', lw=1.5)
    med = float(np.median(x))
    axi.axvline(med, color='red', lw=2, ls='--', label=f'median={med:.3f}')
    axi.set_xlabel(name); axi.set_ylabel("density"); axi.legend()
    axi.set_title(f"marginal: {name}")

ax[1,0].scatter(xy[:, 0], xy[:, 1], s=1, alpha=0.15, c='steelblue')
ax[1,0].scatter([mu1[0], mu2[0]], [mu1[1], mu2[1]], s=120, marker='X',
                c=['#2a6de0', '#d97706'], edgecolor='black', linewidth=1.5,
                label=f'GMM means (w={g2.weights_[0]:.2f}, {g2.weights_[1]:.2f})')
ax[1,0].scatter([joint_median[0]], [joint_median[1]], s=160, marker='+', c='red', linewidth=3,
                label=f'pooled median')
ax[1,0].set_xlabel(p1); ax[1,0].set_ylabel(p2); ax[1,0].legend(fontsize=8)
ax[1,0].set_title(f"joint ({p1}, {p2}) — ΔBIC={joint_bic:.1f}")

t = np.linspace(0, 1, 100)
ax[1,1].plot(t, dens_line, 'k-', lw=2)
ax[1,1].axhline(peaks_min, color='gray', ls=':', label='lower peak height')
ax[1,1].set_xlabel(f"position along line from mean A → mean B")
ax[1,1].set_ylabel("density (GMM k=2)")
ax[1,1].set_title(f"density profile along inter-mode line  (valley/peak = {valley_ratio_2d:.2f})")
ax[1,1].legend(fontsize=8)

plt.tight_layout()
png = os.path.join(out_dir, f"bimodality_{args.chain_prefix}.png")
plt.savefig(png, dpi=110)
print(f"Wrote {png}")
