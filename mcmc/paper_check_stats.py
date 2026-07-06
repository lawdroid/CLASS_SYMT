#!/usr/bin/env python3
"""
paper_check_stats.py — resolve the three v7 CHECK notes in v8_afterglow_dark_energy.tex
that need chain statistics (2026-07-03).

  CHECK @ tex line 1477 (Sec. 5.5, identifiability of c_D):
    Recompute the two-sample Kolmogorov-Smirnov statistic between the 3a and 3b
    c_D marginals on *thinned, approximately independent* draws. The effective
    sample size comes from the integrated autocorrelation time (IAT) of the
    weight-expanded c_D series, computed per chain with the emcee/Sokal FFT
    method. Reports: tau per chain, ESS, D_KS and p-value on the thinned draws,
    plus the old 5e4-resample protocol for continuity with the printed 0.029.

  CHECK @ tex line 1606 (Sec. 5.5, Bayesian model comparison):
    Reconcile the k-NN (Heavens 2017) and Savage-Dickey estimators. Loads the
    2026-06-24 results files (evidence_knn_results.json,
    evidence_vs_lcdm_sddr_results.json) produced against these same chains and
    emits the comparison table. Key structural fact: k-NN yields absolute
    log Z per branch and *inter-branch* deltas only; Delta log Z vs LambdaCDM
    exists only through SDDR (or the staged PolyChord run). The tex paragraph
    must not attribute the {+0.53,+0.35,+0.02} numbers to the k-NN estimator.

  CHECK @ tex line 1196 (Sec. 5.3, cs2 insensitivity run pair):
    Scans the chain tree for any run with cs2_X != 0.5 (e.g. the claimed
    0.25 / 1.20 pair). This script reports what it finds; the box-side hunt is
    driven by run_paper_checks.ps1.

Usage (from CLASS_SYMT root, WSL or the Ubuntu box):
    python3 mcmc/paper_check_stats.py [--repo PATH] [--burn-frac 0.3]
Writes paper_check_results.json next to the repo root.
"""

import argparse, glob, json, os, re, sys
from pathlib import Path

import numpy as np
from scipy.stats import ks_2samp

CHAINS = {
    "3a": ("chains/phase3a_baseline", "afterglow_3a"),
    "3b": ("chains/phase3b_noncrossing", "afterglow_3b"),
}
PARAM = "c_D"  # derived column; KS is invariant under the monotone log10 map anyway


def read_columns(txt, colnames):
    with open(txt) as f:
        header = f.readline().lstrip("#").split()
    idx = [header.index(c) for c in colnames]
    data = np.loadtxt(txt, usecols=idx)
    return data if data.ndim == 2 else data[None, :]


def integrated_autocorr_time(x, c=5.0):
    """Sokal-windowed IAT via FFT autocorrelation (emcee's method)."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n < 100:
        return np.nan
    nfft = 1 << (2 * n - 1).bit_length()
    xc = x - x.mean()
    f = np.fft.rfft(xc, nfft)
    acf = np.fft.irfft(f * np.conjugate(f), nfft)[:n].real
    acf /= acf[0]
    taus = 2.0 * np.cumsum(acf) - 1.0
    window = np.arange(len(taus)) < c * taus
    if window.all():
        return taus[-1]
    m = np.argmin(window)
    return taus[m]


def load_branch(repo, key, burn_frac):
    """Return per-chain weight-expanded c_D series (post burn-in)."""
    d, prefix = CHAINS[key]
    files = sorted(glob.glob(str(repo / d / f"{prefix}.[0-9].txt")))
    if not files:
        raise FileNotFoundError(f"no chain files under {repo / d}")
    expanded = []
    for f in files:
        wc = read_columns(f, ["weight", PARAM])
        nburn = int(burn_frac * len(wc))
        wc = wc[nburn:]
        w = wc[:, 0].astype(int)
        expanded.append(np.repeat(wc[:, 1], w))
    return expanded


def ks_thinned(repo, burn_frac, resample_n=50_000, seed=20260703):
    out = {"param": PARAM, "burn_frac": burn_frac, "branches": {}}
    thinned = {}
    for key in CHAINS:
        chains = load_branch(repo, key, burn_frac)
        taus, draws = [], []
        for x in chains:
            tau = integrated_autocorr_time(x)
            taus.append(tau)
            step = max(1, int(np.ceil(tau)))
            draws.append(x[::step])
        ind = np.concatenate(draws)
        thinned[key] = ind
        n_exp = int(sum(len(x) for x in chains))
        out["branches"][key] = {
            "n_expanded": n_exp,
            "tau_per_chain": [round(float(t), 1) for t in taus],
            "ess_total": int(sum(len(x) / t for x, t in zip(chains, taus))),
            "n_thinned_draws": int(len(ind)),
            "median": float(np.median(ind)),
        }
    D, p = ks_2samp(thinned["3a"], thinned["3b"])
    out["ks_thinned"] = {"D": float(D), "p_value": float(p)}
    out["median_shift_3b_minus_3a"] = round(
        out["branches"]["3b"]["median"] - out["branches"]["3a"]["median"], 3)

    # continuity with the printed protocol: 5e4 equal-weight resampled draws each
    rng = np.random.default_rng(seed)
    res = {}
    for key in CHAINS:
        chains = load_branch(repo, key, burn_frac)
        allx = np.concatenate(chains)
        res[key] = rng.choice(allx, size=resample_n, replace=True)
    D0, p0 = ks_2samp(res["3a"], res["3b"])
    out["ks_resampled_5e4_oldprotocol"] = {"D": float(D0), "nominal_p": float(p0)}
    return out


def evidence_reconciliation(repo):
    knn_f = repo / "evidence_knn_results.json"
    sddr_f = repo / "evidence_vs_lcdm_sddr_results.json"
    if not (knn_f.exists() and sddr_f.exists()):
        return {"error": "evidence results json missing; run evidence_knn.py / evidence_vs_lcdm_sddr.py first"}
    knn = json.loads(knn_f.read_text())
    sddr = json.loads(sddr_f.read_text())
    vs_lcdm = {b: round(sddr[b]["delta_logZ_aft_minus_lcdm"], 3) for b in ("3a", "3b", "3c")}
    sddr_pairs = {
        "3b_minus_3a": vs_lcdm["3b"] - vs_lcdm["3a"],
        "3c_minus_3a": vs_lcdm["3c"] - vs_lcdm["3a"],
        "3c_minus_3b": vs_lcdm["3c"] - vs_lcdm["3b"],
    }
    knn_pairs = {k: round(v["delta_logZ"], 3) for k, v in knn["deltas"].items()}
    agree = {k: round(abs(knn_pairs[k] - round(sddr_pairs[k], 3)), 3) for k in knn_pairs}
    return {
        "sddr_delta_logZ_vs_LCDM": vs_lcdm,
        "sddr_boot_stderr": {b: round(sddr[b]["boot_stderr"], 3) for b in ("3a", "3b", "3c")},
        "knn_absolute_logZ": {b: round(knn["branches"][b]["log_Z"], 3) for b in ("3a", "3b", "3c")},
        "knn_interbranch_deltas": knn_pairs,
        "sddr_interbranch_deltas": {k: round(v, 3) for k, v in sddr_pairs.items()},
        "interbranch_agreement_abs": agree,
        "note": ("k-NN provides absolute log Z per branch and inter-branch deltas; "
                 "Delta log Z vs LCDM comes only from SDDR (LCDM chain not run; "
                 "PolyChord staged). Do not attribute the vs-LCDM numbers to k-NN."),
    }


def cs2_run_hunt(repo):
    hits = []
    for y in glob.glob(str(repo / "chains" / "**" / "*.yaml"), recursive=True) + \
             glob.glob(str(repo / "mcmc" / "*.yaml")):
        try:
            text = Path(y).read_text()
        except OSError:
            continue
        for m in re.finditer(r"cs2_X\s*:\s*([0-9.eE+-]+)", text):
            val = float(m.group(1))
            if abs(val - 0.5) > 1e-9:
                hits.append({"file": os.path.relpath(y, repo), "cs2_X": val})
    return {
        "runs_with_cs2_not_0p5": hits,
        "note": ("empty list = no cs2_X != 0.5 chain/config found in this tree; "
                 "the only recorded sensitivity artifact is scripts/cb2_sensitivity.py "
                 "(spectra-level, cs2 = 1/2 vs 1/3). The 0.25/1.20 pair in the tex "
                 "is unverified unless the box-side hunt finds it."),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=str(Path(__file__).resolve().parent.parent))
    ap.add_argument("--burn-frac", type=float, default=0.3)
    args = ap.parse_args()
    repo = Path(args.repo)

    results = {"generated": "paper_check_stats.py", "repo": str(repo)}
    print("[1/3] KS on ESS-thinned draws (tex line 1477) ...", flush=True)
    results["check_1477_ks_thinned"] = ks_thinned(repo, args.burn_frac)
    print("[2/3] evidence reconciliation (tex line 1606) ...", flush=True)
    results["check_1606_evidence"] = evidence_reconciliation(repo)
    print("[3/3] cs2 run hunt (tex line 1196) ...", flush=True)
    results["check_1196_cs2_runs"] = cs2_run_hunt(repo)

    out = repo / "paper_check_results.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {out}\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
