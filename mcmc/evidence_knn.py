"""
k-nearest-neighbor evidence estimator from Heavens, Sellentin, de Mijolla, Vianello
(2017, arXiv:1704.03472). Reads cobaya chains and returns log Z.

Algorithm:
  - For each weighted sample θ_i in the chain, compute distance d_k(θ_i) to its
    k-th nearest neighbor in the *sampled-parameter* subspace.
  - The k-NN density estimator gives p̂(θ_i) = k / (N_eff · V_d(d_k(i))) where
    V_d(r) is the volume of a d-ball of radius r.
  - The standard importance-sampling identity Z = E_p[L·π/p̂] applied to the
    posterior-distributed samples gives
        Z ≈ (1/k) Σ_i w_i · L_i · π_i · V_d(d_k(i))
        log Z ≈ -log k + logsumexp_i[ -minuslogpost_i + log w_i + d·log d_k(i) + log(π^(d/2)/Γ(d/2+1)) ]
  - We use the cobaya-recorded minuslogpost (= -log(L·π)), so we don't have to
    reconstruct the prior or likelihood separately.

Inputs needed per branch:
  - chain .txt files (cobaya format with header)
  - the .input.yaml so we can identify which columns are *sampled* params
    (only those go into the k-NN distance; derived params would inflate d).

Run from CLASS_SYMT/ root:
    python mcmc/evidence_knn.py
    python mcmc/evidence_knn.py --k 10 --burn-frac 0.4

Caveats: k-NN evidence is known to bias high in high dimensions; Δ log Z
between two chains with the same sampler tend to cancel out systematics, so
the cross-check Δ log Z(3b − 3a) ≈ 0 (predicted from flat-β analysis) is a
more reliable diagnostic than the absolute values.
"""

import argparse, json, math, os, re, sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

BRANCHES = {
    "3a": "chains/phase3a_baseline",
    "3b": "chains/phase3b_noncrossing",
    "3c": "chains/phase3c_b1narrow",
}


def sampled_params_from_yaml(yaml_path: Path):
    """Identify sampled params: in cobaya, sampled => has prior:{min,max}; derived
    or fixed => doesn't. We do a lightweight regex pass to avoid pulling a yaml
    library dependency."""
    text = yaml_path.read_text()
    sampled = []
    in_params = False
    cur_param = None
    cur_block = []
    indent_param = None

    for line in text.splitlines():
        if re.match(r"^params:\s*$", line):
            in_params = True
            continue
        if not in_params:
            continue
        if re.match(r"^\S", line) and not line.startswith("params:"):
            in_params = False
            continue
        m = re.match(r"^  ([A-Za-z_][\w]*):\s*$", line)
        if m:
            if cur_param is not None:
                if any("prior:" in l for l in cur_block) and not any("value:" in l or "derived:" in l for l in cur_block):
                    sampled.append(cur_param)
            cur_param = m.group(1); cur_block = []
            continue
        if cur_param is not None:
            cur_block.append(line)
    if cur_param is not None:
        if any("prior:" in l for l in cur_block) and not any("value:" in l or "derived:" in l for l in cur_block):
            sampled.append(cur_param)
    return sampled


def read_chain(chain_dir: Path, glob_pat: str, sampled_cols, burn_frac):
    """Return weights array (n,), minuslogpost array (n,), points array (n, d)."""
    import numpy as np
    files = sorted(chain_dir.glob(glob_pat))
    if not files:
        return None
    with files[0].open() as fh:
        header = fh.readline().lstrip("#").strip().split()
    weights, mlpost, points = [], [], []
    col_idx = {c: i for i, c in enumerate(header)}
    missing = [c for c in sampled_cols if c not in col_idx]
    if missing:
        print(f"  WARNING: sampled-param columns not in chain header: {missing}")
        sampled_cols = [c for c in sampled_cols if c in col_idx]
    sp_idx = [col_idx[c] for c in sampled_cols]
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
                weights.append(float(vals[0]))
                mlpost.append(float(vals[1]))
                points.append([float(vals[i]) for i in sp_idx])
            except ValueError:
                continue
    if not weights:
        return None
    return (np.asarray(weights), np.asarray(mlpost),
            np.asarray(points), sampled_cols)


def log_volume_ball(d, log_r):
    """log volume of d-ball of radius r: d log r + (d/2) log π - log Γ(d/2+1)"""
    return d * log_r + (d/2.0) * math.log(math.pi) - math.lgamma(d/2.0 + 1)


def heavens_logZ(weights, mlpost, points, k=5):
    """k-NN evidence estimator. Returns (log_Z, stderr, n_eff, d)."""
    import numpy as np
    from sklearn.neighbors import NearestNeighbors
    n, d = points.shape

    # Whiten the parameter space so the k-NN distance is sensible across
    # heterogeneously-scaled params (H0 ~ O(70), omega_b ~ O(0.02), etc.).
    mu  = (points * weights[:, None]).sum(0) / weights.sum()
    var = (weights[:, None] * (points - mu)**2).sum(0) / weights.sum()
    std = np.sqrt(var)
    std[std < 1e-12] = 1.0
    X = (points - mu) / std

    # k-NN distances; the chain is weighted, but sklearn doesn't natively
    # weight density. Use unweighted k-NN; combine with weights only in the
    # final estimator sum. (Standard for Cobaya's near-unity weights.)
    nn = NearestNeighbors(n_neighbors=k+1).fit(X)
    dists, _ = nn.kneighbors(X)
    dk = dists[:, k]  # k-th NN distance (excluding self at index 0)
    dk = np.maximum(dk, 1e-300)
    log_dk = np.log(dk)

    # Whitening rescales volumes; need to correct back to original metric:
    # log V_original = log V_whitened + Σ log std_j  (since vol scales as Πstd_j).
    log_vol_correction = float(np.log(std).sum())
    log_vol_d = np.array([log_volume_ball(d, lr) for lr in log_dk]) + log_vol_correction

    # Z = (1/k) Σ w_i L_i π_i V_d(d_k(i)),  log L_i π_i = -mlpost_i
    log_w = np.log(np.maximum(weights, 1e-300))
    terms = -mlpost + log_w + log_vol_d
    # log-sum-exp for stability
    m = terms.max()
    logZ = m + math.log(np.exp(terms - m).sum()) - math.log(k)

    # Crude stderr estimate: split-half
    half = n // 2
    def logZ_subset(idx):
        sub_terms = -mlpost[idx] + np.log(np.maximum(weights[idx], 1e-300)) + log_vol_d[idx]
        m_ = sub_terms.max()
        return m_ + math.log(np.exp(sub_terms - m_).sum()) - math.log(k)
    lz1 = logZ_subset(slice(0, half))
    lz2 = logZ_subset(slice(half, n))
    stderr = abs(lz1 - lz2) / 2.0   # rough — not a real bootstrap

    n_eff = (weights.sum()**2) / (weights**2).sum()
    return logZ, stderr, n_eff, d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--burn-frac", type=float, default=0.3)
    ap.add_argument("--out", default=str(REPO / "evidence_knn_results.json"))
    args = ap.parse_args()

    results = {
        "method": "Heavens 2017 k-NN evidence (arXiv:1704.03472)",
        "k":      args.k,
        "burn_frac": args.burn_frac,
        "branches": {},
    }

    for br, cdir in BRANCHES.items():
        chain_dir = REPO / cdir
        yaml_files = list(chain_dir.glob("*.input.yaml"))
        if not chain_dir.exists() or not yaml_files:
            print(f"== {br}: no chain dir / yaml at {chain_dir}; skip")
            continue
        sampled = sampled_params_from_yaml(yaml_files[0])
        print(f"== {br}: sampled params (d={len(sampled)}): {sampled}")
        glob_pat = f"afterglow_{br}.[1-9].txt"
        chain = read_chain(chain_dir, glob_pat, sampled, args.burn_frac)
        if chain is None:
            print(f"   no chain rows; skip")
            continue
        weights, mlpost, points, sampled_actual = chain
        print(f"   n_samples={len(weights)}, n_eff={(weights.sum()**2)/(weights**2).sum():.0f}, d_effective={points.shape[1]}")
        logZ, stderr, n_eff, d = heavens_logZ(weights, mlpost, points, k=args.k)
        print(f"   log Z = {logZ:.3f} ± {stderr:.3f} (split-half)")
        results["branches"][br] = {
            "log_Z":       logZ,
            "stderr":      stderr,
            "n_eff":       n_eff,
            "d_sampled":   d,
            "sampled_params": sampled_actual,
        }

    # Deltas (consistency cross-checks)
    branches = list(results["branches"].keys())
    deltas = {}
    for i in range(len(branches)):
        for j in range(i+1, len(branches)):
            a, b = branches[i], branches[j]
            la = results["branches"][a]["log_Z"]
            lb = results["branches"][b]["log_Z"]
            sa = results["branches"][a]["stderr"]
            sb = results["branches"][b]["stderr"]
            deltas[f"{b}_minus_{a}"] = {
                "delta_logZ": lb - la,
                "stderr":     math.sqrt(sa*sa + sb*sb),
            }
    results["deltas"] = deltas

    with open(args.out, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nwrote {args.out}")
    for k, v in deltas.items():
        print(f"  Δ log Z ({k}) = {v['delta_logZ']:+.3f} ± {v['stderr']:.3f}")


if __name__ == "__main__":
    main()
