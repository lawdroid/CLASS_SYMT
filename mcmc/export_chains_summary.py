"""
Real-chain export for phase3_summary.html / phase3_branches_explorer.html.
Replaces gen_chains_summary_mock.py once real chains have been fetched.

Reads cobaya chain .txt files from:
  chains/phase3a_baseline/afterglow_3a.[1-4].txt
  chains/phase3b_noncrossing/afterglow_3b.[1-4].txt
  chains/phase3c_b1narrow/afterglow_3c.[1-4].txt

Writes chains_summary.json with:
  - data_source: "REAL chains, <timestamp>"
  - per-branch stats_1d, kde_2d (same schema as mock)
  - pending_numbers: {beta_lower95_3b, logZ_3a, logZ_3b, delta_logZ_3b_vs_lcdm}

logZ values are computed via MCEvidence if importable; otherwise null.
delta_logZ requires a ΛCDM chain at chains/lcdm_reference/; null otherwise.

Run from CLASS_SYMT/ root:
  python mcmc/export_chains_summary.py
  python mcmc/export_chains_summary.py --burn-frac 0.4 --kde-grid 60
"""

import argparse, glob, json, math, os, sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

BRANCH_CONFIG = {
    "3a": {
        "chain_dir":  "chains/phase3a_baseline",
        "chain_glob": "afterglow_3a.[1-9].txt",
        "label": "3a baseline (β = 0, arrested-medium)",
        "color": "#155fa5",
        "beta_fixed": 0.0,
    },
    "3b": {
        "chain_dir":  "chains/phase3b_noncrossing",
        "chain_glob": "afterglow_3b.[1-9].txt",
        "label": "3b non-crossing (β ∈ (0, 2.4])",
        "color": "#0f6e56",
        "beta_fixed": None,
    },
    "3c": {
        "chain_dir":  "chains/phase3c_b1narrow",
        "chain_glob": "afterglow_3c.[1-9].txt",
        "label": "3c crossing, B1-narrow (β ∈ [2.65, 2.75])",
        "color": "#993c1d",
        "beta_fixed": None,
    },
}

# Chain column → exported param name
COL_MAP = {
    "beta_aft": "beta",
    "H0":       "H0",
    "Omega_m":  "Omega_m",
    "sigma8":   "sigma8",
    "w_X":      "w_chi",
    "c_D":      "c_D",
}
# h is H0/100 — derived in-script for compatibility with the mock schema
PARAMS = ["beta", "Omega_m", "h", "sigma8", "S8", "w_chi", "c_D"]
PRETTY = {"beta":"β","Omega_m":"Ωₘ","h":"h","sigma8":"σ₈","S8":"S₈","w_chi":"w_χ","c_D":"c_D"}
PAIRS = [
    ("beta", "c_D"), ("beta", "Omega_m"), ("beta", "sigma8"),
    ("beta", "h"), ("beta", "w_chi"),
    ("Omega_m", "sigma8"), ("Omega_m", "h"), ("h", "sigma8"),
    ("Omega_m", "S8"), ("c_D", "sigma8"), ("c_D", "Omega_m"),
    ("w_chi", "c_D"),
]


def read_cobaya_chain(chain_dir: Path, glob_pat: str, burn_frac: float):
    """Return (weights, samples_dict, minus_log_post_list) or None if no files."""
    files = sorted(chain_dir.glob(glob_pat))
    if not files:
        return None
    with files[0].open() as fh:
        header = fh.readline().lstrip("#").strip().split()
    weights, mlogpost = [], []
    cols = {c: [] for c in header}
    for fp in files:
        with fp.open() as fh:
            lines = fh.readlines()
        rows = [l for l in lines if l.strip() and not l.lstrip().startswith("#")]
        burn = int(len(rows) * burn_frac)
        for line in rows[burn:]:
            vals = line.split()
            if len(vals) != len(header):
                continue
            try:
                fvals = [float(v) for v in vals]
            except ValueError:
                continue
            for c, v in zip(header, fvals):
                cols[c].append(v)
            weights.append(fvals[0])
            mlogpost.append(fvals[1])
    if not weights:
        return None
    return weights, cols, mlogpost


def build_samples(branch_key, chain_data):
    """Project chain columns onto the export schema (PARAMS)."""
    weights, cols, _ = chain_data
    n = len(weights)
    samples = []
    beta_fixed = BRANCH_CONFIG[branch_key]["beta_fixed"]
    for i in range(n):
        rec = {}
        if beta_fixed is not None:
            rec["beta"] = beta_fixed
        else:
            rec["beta"] = cols.get("beta_aft", [0.0])[i]
        rec["Omega_m"] = cols["Omega_m"][i]
        rec["sigma8"]  = cols["sigma8"][i]
        rec["w_chi"]   = cols["w_X"][i]
        rec["c_D"]     = cols["c_D"][i]
        rec["h"]       = cols["H0"][i] / 100.0
        rec["S8"]      = rec["sigma8"] * math.sqrt(rec["Omega_m"] / 0.3)
        samples.append(rec)
    return samples, weights


def w_quantile(xs, ws, q):
    pairs = sorted(zip(xs, ws), key=lambda t: t[0])
    total = sum(ws)
    cum = 0.0
    for x, w in pairs:
        cum += w
        if cum >= q * total:
            return x
    return pairs[-1][0]


def w_mean(xs, ws):
    return sum(x*w for x, w in zip(xs, ws)) / sum(ws)


def w_std(xs, ws):
    m = w_mean(xs, ws)
    return math.sqrt(sum(w*(x-m)**2 for x, w in zip(xs, ws)) / sum(ws))


def stats_1d(samples, weights, key):
    xs = [s[key] for s in samples]
    return {
        "median": w_quantile(xs, weights, 0.5),
        "q16":    w_quantile(xs, weights, 0.16),
        "q84":    w_quantile(xs, weights, 0.84),
        "mean":   w_mean(xs, weights),
        "std":    w_std(xs, weights),
        "min":    min(xs),
        "max":    max(xs),
    }


def kde2d(samples, weights, xk, yk, nx=50, ny=50):
    xs = [s[xk] for s in samples]
    ys = [s[yk] for s in samples]
    xmn, xmx = min(xs), max(xs)
    ymn, ymx = min(ys), max(ys)
    if xmx - xmn < 1e-4:
        c = (xmx+xmn)/2; xmn, xmx = c-0.04, c+0.04
    if ymx - ymn < 1e-4:
        c = (ymx+ymn)/2; ymn, ymx = c-0.04, c+0.04
    padx = 0.06*(xmx-xmn); pady = 0.06*(ymx-ymn)
    xmin = xmn - padx; xmax = xmx + padx
    ymin = ymn - pady; ymax = ymx + pady
    sx = max(w_std(xs, weights), (xmax-xmin)*0.04, 1e-4)
    sy = max(w_std(ys, weights), (ymax-ymin)*0.04, 1e-4)
    n_eff = sum(weights)**2 / sum(w*w for w in weights)
    hx = 1.06 * sx * (n_eff ** -0.2)
    hy = 1.06 * sy * (n_eff ** -0.2)
    dx = (xmax-xmin)/(nx-1); dy = (ymax-ymin)/(ny-1)
    grid = [[0.0]*nx for _ in range(ny)]
    rx = max(2, int(math.ceil(3*hx/dx)))
    ry = max(2, int(math.ceil(3*hy/dy)))
    for i_s in range(len(xs)):
        sxv = xs[i_s]; syv = ys[i_s]; w = weights[i_s]
        ix = int(round((sxv-xmin)/dx)); iy = int(round((syv-ymin)/dy))
        i0 = max(0, ix-rx); i1 = min(nx-1, ix+rx)
        j0 = max(0, iy-ry); j1 = min(ny-1, iy+ry)
        for j in range(j0, j1+1):
            gy = ymin + j*dy
            ty = (gy-syv)/hy
            ky = math.exp(-0.5*ty*ty)
            if ky < 0.002: continue
            for i in range(i0, i1+1):
                gx = xmin + i*dx
                tx = (gx-sxv)/hx
                grid[j][i] += w * math.exp(-0.5*tx*tx) * ky
    def r4(v): return 0 if v == 0 else float(f"{v:.4g}")
    flat = [r4(grid[j][i]) for j in range(ny) for i in range(nx)]
    return {
        "xmin": round(xmin, 5), "xmax": round(xmax, 5),
        "ymin": round(ymin, 5), "ymax": round(ymax, 5),
        "nx": nx, "ny": ny, "grid": flat,
    }


def gelman_rubin(chain_dir: Path, glob_pat: str, key="c_D", burn_frac=0.3):
    """Crude per-chain R-1 on a single param; used as a sanity badge for the page."""
    files = sorted(chain_dir.glob(glob_pat))
    if len(files) < 2: return None
    means, vars_, ns = [], [], []
    for fp in files:
        with fp.open() as fh:
            header = fh.readline().lstrip("#").strip().split()
            if key not in header: return None
            idx = header.index(key)
            wi = 0
            xs, ws = [], []
            for line in fh:
                if not line.strip() or line.startswith("#"): continue
                vals = line.split()
                if len(vals) != len(header): continue
                ws.append(float(vals[0])); xs.append(float(vals[idx]))
        burn = int(len(xs)*burn_frac)
        xs = xs[burn:]; ws = ws[burn:]
        if not xs: return None
        m = w_mean(xs, ws); v = w_std(xs, ws)**2
        means.append(m); vars_.append(v); ns.append(sum(ws))
    n_avg = sum(ns)/len(ns)
    W = sum(vars_)/len(vars_)
    mean_of_means = sum(means)/len(means)
    B = n_avg * sum((m-mean_of_means)**2 for m in means) / (len(means)-1)
    if W <= 0: return None
    Vhat = (1 - 1/n_avg)*W + B/n_avg
    return math.sqrt(Vhat/W) - 1.0


def try_mcevidence(chain_dir: Path, glob_pat: str, burn_frac=0.3):
    """Return log Z estimate via MCEvidence if available; None otherwise."""
    try:
        from MCEvidence import MCEvidence  # heavens-et-al-2017
    except ImportError:
        try:
            import mcevidence as _mc
            MCEvidence = _mc.MCEvidence
        except Exception:
            return None
    try:
        files = sorted(chain_dir.glob(glob_pat))
        if not files: return None
        root = str(files[0]).rsplit(".", 2)[0]  # strip ".N.txt"
        mce = MCEvidence(root, burnlen=burn_frac, ndim=None)
        logZ = mce.evidence()
        if isinstance(logZ, (list, tuple)):
            logZ = float(logZ[0])
        return float(logZ)
    except Exception as e:
        print(f"  MCEvidence failed: {e}", file=sys.stderr)
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--burn-frac", type=float, default=0.3)
    ap.add_argument("--kde-grid",  type=int,   default=50)
    ap.add_argument("--out", default=str(REPO/"chains_summary.json"))
    args = ap.parse_args()

    import datetime
    out = {
        "schema_version": 1,
        "generated_utc":  datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "data_source":    "REAL chains",
        "params":         PARAMS,
        "pretty":         PRETTY,
        "credible_levels":[0.68, 0.95],
        "branches":       {},
        "pending_numbers":{
            "beta_lower95_3b":       None,
            "logZ_3a":               None,
            "logZ_3b":               None,
            "delta_logZ_3b_vs_lcdm": None,
        },
    }

    missing_branches = []
    for br, cfg in BRANCH_CONFIG.items():
        chain_dir = REPO / cfg["chain_dir"]
        print(f"== {br}: {chain_dir}")
        chain = read_cobaya_chain(chain_dir, cfg["chain_glob"], args.burn_frac)
        if chain is None:
            print(f"   NO CHAIN FILES — skipping {br}")
            missing_branches.append(br)
            continue
        samples, weights = build_samples(br, chain)
        print(f"   {len(samples)} weighted samples")

        R1 = gelman_rubin(chain_dir, cfg["chain_glob"], "c_D", args.burn_frac)
        branch_data = {
            "label":     cfg["label"],
            "color":     cfg["color"],
            "R1":        R1,
            "n_samples": len(samples),
            "stats_1d":  {p: stats_1d(samples, weights, p) for p in PARAMS},
            "kde_2d":    {},
        }
        for xk, yk in PAIRS:
            branch_data["kde_2d"][f"{xk}__{yk}"] = kde2d(
                samples, weights, xk, yk, nx=args.kde_grid, ny=args.kde_grid)
        out["branches"][br] = branch_data

        if br == "3b":
            beta_xs = [s["beta"] for s in samples]
            out["pending_numbers"]["beta_lower95_3b"] = w_quantile(beta_xs, weights, 0.025)

        if br in ("3a", "3b"):
            logZ = try_mcevidence(chain_dir, cfg["chain_glob"], args.burn_frac)
            if logZ is not None:
                out["pending_numbers"][f"logZ_{br}"] = logZ
                print(f"   log Z ({br}) = {logZ:.3f}")

    # delta log Z requires a ΛCDM reference chain
    lcdm_dir = REPO / "chains/lcdm_reference"
    if lcdm_dir.exists():
        logZ_lcdm = try_mcevidence(lcdm_dir, "*.[1-9].txt", args.burn_frac)
        lz3b = out["pending_numbers"]["logZ_3b"]
        if logZ_lcdm is not None and lz3b is not None:
            out["pending_numbers"]["delta_logZ_3b_vs_lcdm"] = lz3b - logZ_lcdm
            out["pending_numbers"]["logZ_lcdm"] = logZ_lcdm

    if missing_branches:
        out["data_source"] += f" (branches missing: {','.join(missing_branches)})"

    # Safety: refuse to clobber an existing chains_summary.json with an empty
    # export when no real chains were found. This protects the mock during
    # dry runs before fetch_phase3_chains_for_tom.ps1 has been executed.
    if len(out["branches"]) == 0 and os.path.exists(args.out):
        print("\nNo real chains found and chains_summary.json already exists — "
              "refusing to overwrite. Run fetch_phase3_chains_for_tom.ps1 first.")
        sys.exit(1)

    with open(args.out, "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    size_kib = os.path.getsize(args.out) / 1024
    print(f"\nwrote {args.out}  ({size_kib:.1f} KiB)")
    print(f"data_source: {out['data_source']}")
    print("pending_numbers:")
    for k, v in out["pending_numbers"].items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
