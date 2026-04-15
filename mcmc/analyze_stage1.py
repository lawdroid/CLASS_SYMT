#!/usr/bin/env python3
"""
analyze_stage1.py — post-Stage-1 triage.

Loads MCMC chains, prints Gelman-Rubin, marginals, correlations with
c_D and beta_aft, and writes a triangle plot.

Usage:  python analyze_stage1.py
"""
import os
import sys
from pathlib import Path

try:
    from getdist.mcsamples import loadMCSamples
    from getdist import plots
    import numpy as np
except ImportError as e:
    print(f"[ERROR] Missing package: {e}. Install: pip install getdist numpy matplotlib")
    sys.exit(1)

CHAINS = Path(__file__).parent / "chains" / "stage1"
ROOT = CHAINS / "afterglow_bg"

if not any(CHAINS.glob("afterglow_bg*.txt")):
    print(f"[ERROR] No chains at {ROOT}*.txt")
    print("        Run submit_stage1.{ps1,sh} first.")
    sys.exit(1)

s = loadMCSamples(str(ROOT))

# ----- convergence -----
print("=" * 60)
print("  Stage 1 — Gelman-Rubin R-1 diagnostic")
print("=" * 60)
print(s.getConvergeTests(what=['MeanVar', 'CorrLengths', 'GelmanRubin']))

# ----- marginals -----
print("\n" + "=" * 60)
print("  1D marginals (mean ± 68% CL)")
print("=" * 60)
for p in ["H0", "omega_b", "omega_cdm", "c_D", "beta_aft", "Omega_X", "w_X"]:
    try:
        stats = s.getMargeStats().parWithName(p)
        print(f"  {p:10s} = {stats.mean:+.4f}  ±  {stats.err:.4f}   (68% CL)")
    except Exception:
        pass

# ----- headline afterglow result -----
print("\n" + "=" * 60)
print("  Afterglow summary")
print("=" * 60)
try:
    cD  = s.getMargeStats().parWithName("c_D")
    bA  = s.getMargeStats().parWithName("beta_aft")
    wX  = s.getMargeStats().parWithName("w_X")
    print(f"  c_D     = {cD.mean:.3f} ± {cD.err:.3f}")
    print(f"  beta    = {bA.mean:.3f} ± {bA.err:.3f}")
    print(f"  w_X(0)  = {wX.mean:.3f} ± {wX.err:.3f}")
    if wX.mean + wX.err < -0.97:
        print("  => afterglow indistinguishable from LCDM at 1 sigma")
    elif wX.mean + 2 * wX.err < -1.0:
        print("  => w_X > -1 detection: DYNAMICAL DARK ENERGY at >2 sigma")
except Exception as e:
    print(f"  (summary failed: {e})")

# ----- triangle plot -----
print("\n[INFO] Writing triangle plot to chains/stage1/triangle.pdf ...")
g = plots.get_subplot_plotter(width_inch=10)
g.triangle_plot(s, ["H0", "omega_cdm", "c_D", "beta_aft", "w_X"],
                filled=True)
g.export(str(CHAINS / "triangle.pdf"))
print("       done.")
