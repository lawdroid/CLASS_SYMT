"""
Delta log Z (afterglow branch - LCDM) via the Savage-Dickey density ratio at
the nested LCDM limit.

LCDM is the c_D -> infinity limit of every afterglow branch:
    w_X = -1 + 1/(3 c_D)  ->  -1,   and the interaction (beta/c_D)*Psi -> 0.
So LCDM sits at the finite point  u == 1 + w_X = 1/(3 c_D) = 0.

For a nested model, the Bayes factor is the SDDR
    B_{LCDM,aft} = Z_LCDM / Z_aft = p(u=0 | d) / pi(u=0)
where p(.|d) and pi(.) are the *marginal* posterior and prior densities in u.
(beta marginalizes out: at u=0 the interaction vanishes and beta is unidentified,
so the same u-marginal SDDR applies to 3a, 3b, 3c.)

The ratio r(u) = p(u|d)/pi(u) is the normalized profile likelihood, finite at
u->0 even though pi(u) ~ 1/u diverges (the divergences cancel). We estimate
r(u) on the prior support and extrapolate r(0); several extrapolation forms
bracket the systematic.

    Delta log Z (aft - LCDM) = -log B_{LCDM,aft} = -log r(0).

Induced prior on u from the flat prior on L=log10_c_D over [0.7,2.0]:
    u = 1/(3*10^L)  =>  pi(u) = (1/1.3) * 1/(u ln10),  u in [1/300, 1/15].
"""
import json, glob, numpy as np
from scipy.stats import gaussian_kde

LN10 = np.log(10.0)
L_LO, L_HI = 0.7, 2.0           # flat prior on log10_c_D
U_LO, U_HI = 1/(3*10**L_HI), 1/(3*10**L_LO)   # u range [0.003333, 0.066667]
BURN = 0.3

BRANCHES = {
    "3a": ("chains/phase3a_baseline/afterglow_3a.[1-4].txt", 21),   # w_X col
    "3b": ("chains/phase3b_noncrossing/afterglow_3b.[1-4].txt", 22),
    "3c": ("chains/phase3c_b1narrow/afterglow_3c.[1-4].txt", 22),
}

def load_u(glob_pat, wx_col):
    W, U = [], []
    for f in sorted(glob.glob(glob_pat)):
        d = np.loadtxt(f, usecols=(0, wx_col))
        n = len(d); d = d[int(BURN*n):]
        W.append(d[:, 0]); U.append(1.0 + d[:, 1])   # u = 1 + w_X
    return np.concatenate(W), np.concatenate(U)

def prior_u(u):
    return (1.0/(L_HI-L_LO)) / (u*LN10)

def sddr_logZ(W, U):
    # weighted KDE of posterior in u; restrict to prior support
    m = (U >= U_LO-1e-9) & (U <= U_HI+1e-9)
    W, U = W[m], U[m]
    kde = gaussian_kde(U, weights=W)
    # grid over support; r(u) = p(u|d)/pi(u) = normalized profile likelihood
    g = np.linspace(U_LO, U_HI, 600)
    r = kde(g) / prior_u(g)
    # renormalize so prior-average of r is 1 (E_pi[r]=1) -- guards KDE leakage
    pi = prior_u(g); dz = g[1]-g[0]
    norm = np.sum(r*pi)*dz
    r = r/norm
    # extrapolate r(u)->u=0 with several forms over the lower third
    lo = g <= (U_LO + (U_HI-U_LO)/3)
    x, y = g[lo], r[lo]
    fits = {}
    # linear
    p1 = np.polyfit(x, y, 1); fits["linear"]   = np.polyval(p1, 0.0)
    # quadratic
    p2 = np.polyfit(x, y, 2); fits["quadratic"] = np.polyval(p2, 0.0)
    # log-linear (exp fit) guards positivity
    pl = np.polyfit(x, np.log(np.maximum(y,1e-12)), 1); fits["loglinear"] = np.exp(np.polyval(pl, 0.0))
    # flat (value at lowest edge) -- conservative floor
    fits["edge"] = float(y[0])
    r0_vals = np.array([fits["linear"], fits["quadratic"], fits["loglinear"], fits["edge"]])
    r0 = np.median(r0_vals)
    dlogZ = -np.log(r0)                  # Delta log Z (aft - LCDM)
    spread = np.array([-np.log(max(v,1e-9)) for v in r0_vals])
    return dlogZ, fits, r0, (spread.min(), spread.max()), len(U)

# bootstrap error on r0 (resample chains' weighted samples)
def sddr_boot(W, U, nboot=40, seed=12345):
    rng = np.random.default_rng(seed)
    m = (U >= U_LO-1e-9) & (U <= U_HI+1e-9)
    W, U = W[m], U[m]
    p = W/ W.sum()
    n = min(len(U), 60000)               # subsample for speed
    out = []
    for _ in range(nboot):
        idx = rng.choice(len(U), size=n, replace=True, p=p)
        u = U[idx]
        kde = gaussian_kde(u)
        g = np.linspace(U_LO, U_HI, 400); pi = prior_u(g)
        r = kde(g)/pi; r/= np.sum(r*pi)*(g[1]-g[0])
        lo = g <= (U_LO + (U_HI-U_LO)/3)
        r0 = np.polyval(np.polyfit(g[lo], r[lo], 2), 0.0)
        if r0>0: out.append(-np.log(r0))
    return np.std(out)

results = {}
for name,(pat,wcol) in BRANCHES.items():
    W,U = load_u(pat, wcol)
    dlogZ, fits, r0, (smin,smax), nu = sddr_logZ(W,U)
    berr = sddr_boot(W,U)
    results[name] = dict(delta_logZ_aft_minus_lcdm=dlogZ, r0=r0,
                         extrap_min=smin, extrap_max=smax, boot_stderr=berr,
                         fits={k:float(v) for k,v in fits.items()},
                         u_median=float(np.median(U)), n=nu)
    print(f"{name}: dlogZ(aft-LCDM)={dlogZ:+.3f}  (extrap [{smin:+.3f},{smax:+.3f}], boot±{berr:.3f})  r0={r0:.3f}")

# cross-check vs kNN inter-branch deltas
kg = json.load(open("evidence_knn_results.json"))["deltas"]
print("\n=== cross-check: SDDR inter-branch differences vs kNN ===")
def d(a,b): return results[a]["delta_logZ_aft_minus_lcdm"]-results[b]["delta_logZ_aft_minus_lcdm"]
for pair,key in [(("3b","3a"),"3b_minus_3a"),(("3c","3a"),"3c_minus_3a"),(("3c","3b"),"3c_minus_3b")]:
    sd=d(*pair); kn=kg[key]["delta_logZ"]
    print(f"  {pair[0]}-{pair[1]}: SDDR {sd:+.3f}   kNN {kn:+.3f}   diff {sd-kn:+.3f}")

json.dump(results, open("evidence_vs_lcdm_sddr_results.json","w"), indent=2)
print("\nsaved evidence_vs_lcdm_sddr_results.json")
