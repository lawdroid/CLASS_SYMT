#!/usr/bin/env python3
"""
forward_model_sweep.py — Phase 3 pre-flight beta sweep.

Per Tom/Phase3_MCMC_Windows_Execution.md §5.1: runs CLASS_SYMT at
seven beta values spanning the three-branch design (3a baseline,
3b non-crossing strip, 3c crossing strip), at fixed c_D = 8.0,
cs2_X = 0.5.

  beta = 0.0     -> 3a baseline
  beta = 0.5     -> 3b lower
  beta = 1.0     -> 3b mid
  beta = 2.0     -> 3b upper-mid
  beta = 2.7     -> 3c lower
  beta = 3.0     -> 3c mid
  beta = 3.4     -> 3c upper

Plots:
  Top row  : Delta C_l / C_l^LCDM vs l for TT / EE
  Bottom   : log P(k) absolute curves overlaid (LCDM + 7 betas)
  Saved as figures/forward_model_sweep.pdf

Numerical summary printed to stdout: max |Delta C_l / C_l| per beta
across the masked l range, plus relative P(k) at k = 0.1 h/Mpc.

Usage (from CLASS_SYMT root):
    python3 scripts/forward_model_sweep.py
"""
from __future__ import annotations
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- configuration ----------
CLASS_ROOT = Path(__file__).resolve().parent.parent
CLASS_BIN = CLASS_ROOT / "class"
WORK_DIR = Path("/tmp/forward_model_sweep_work")
FIGURES_DIR = CLASS_ROOT / "figures"
OUTPUT_PDF = FIGURES_DIR / "forward_model_sweep.pdf"

BETAS = [0.0, 0.5, 1.0, 2.0, 2.7, 3.0, 3.4]
C_D = 8.0
CS2_X = 0.5

BASE_INI = """
h = 0.6732
omega_b = 0.02238
omega_cdm = 0.1201
A_s = 2.1e-9
n_s = 0.966
tau_reio = 0.0543
output = tCl,pCl,lCl,mPk
lensing = yes
l_max_scalars = 2500
P_k_max_h/Mpc = 1.0
"""


def write_ini(case_idx: int, beta: float | None) -> Path:
    root = WORK_DIR / f"case_{case_idx:02d}_"
    lines = [BASE_INI, f"root = {root}"]
    if beta is not None:
        lines += [
            "afterglow_on = 1",
            f"c_D = {C_D}",
            f"beta_aft = {beta}",
            f"cs2_X = {CS2_X}",
        ]
    ini_path = WORK_DIR / f"case_{case_idx:02d}.ini"
    ini_path.write_text("\n".join(lines) + "\n")
    return ini_path


def run_class(ini_path: Path) -> None:
    print(f"  running CLASS on {ini_path.name} ...", flush=True)
    proc = subprocess.run(
        [str(CLASS_BIN), str(ini_path)],
        cwd=CLASS_ROOT,
        capture_output=True,
        text=True,
        timeout=300,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise RuntimeError(f"CLASS failed on {ini_path}")


def read_cl_lensed(case_idx: int) -> dict:
    arr = np.loadtxt(WORK_DIR / f"case_{case_idx:02d}_00_cl_lensed.dat")
    return {"l": arr[:, 0], "TT": arr[:, 1], "EE": arr[:, 2], "TE": arr[:, 3]}


def read_pk(case_idx: int) -> dict:
    arr = np.loadtxt(WORK_DIR / f"case_{case_idx:02d}_00_pk.dat")
    return {"k": arr[:, 0], "P": arr[:, 1]}


def main() -> int:
    if not CLASS_BIN.exists():
        sys.stderr.write(f"ERROR: CLASS binary not found at {CLASS_BIN}.\n")
        return 1

    if WORK_DIR.exists():
        shutil.rmtree(WORK_DIR)
    WORK_DIR.mkdir(parents=True)
    FIGURES_DIR.mkdir(exist_ok=True)

    # Case 0 = LCDM reference, cases 1..7 = the seven betas
    cases = [{"label": "LCDM", "beta": None}] + [
        {"label": rf"$\beta={b}$", "beta": b} for b in BETAS
    ]

    print(f"Running {len(cases)} CLASS configurations:")
    for i, c in enumerate(cases):
        print(f"  case {i}: {c['label']}")
        run_class(write_ini(i, c["beta"]))

    cls_data = [read_cl_lensed(i) for i in range(len(cases))]
    pks = [read_pk(i) for i in range(len(cases))]
    lcdm_cl = cls_data[0]
    lcdm_pk = pks[0]

    # Detect cases with non-finite output (CLASS perturbation integrator
    # can fail in the 3c crossing strip at c_D = 8.0). Track but still
    # plot what's plottable.
    bad_cases = []
    for i in range(1, len(cases)):
        cl_bad = not (np.isfinite(cls_data[i]["TT"]).all() and np.isfinite(cls_data[i]["EE"]).all())
        pk_bad = not np.isfinite(pks[i]["P"]).all()
        if cl_bad or pk_bad:
            bad_cases.append((BETAS[i - 1], "Cl" if cl_bad else "", "Pk" if pk_bad else ""))
    if bad_cases:
        print("\nWARNING: numerical issues detected (NaN/Inf in CLASS output):")
        for beta, cl, pk in bad_cases:
            print(f"  beta={beta}: {' '.join(x for x in (cl, pk) if x)}")
        print("Affected curves are skipped in the plot but reported in the summary.\n")

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    cmap = plt.get_cmap("viridis")
    colors = [cmap(x) for x in np.linspace(0.05, 0.95, len(BETAS))]

    def finite_only(*arrs):
        m = np.ones_like(arrs[0], dtype=bool)
        for a in arrs:
            m &= np.isfinite(a)
        return m

    # --- top-left: Delta C^TT / C^TT ---
    ax = axes[0, 0]
    mask_tt = np.abs(lcdm_cl["TT"]) > 0.05 * np.abs(lcdm_cl["TT"]).max()
    for i, b in enumerate(BETAS):
        ratio = (cls_data[i + 1]["TT"] - lcdm_cl["TT"]) / lcdm_cl["TT"]
        m = mask_tt & finite_only(ratio)
        if m.sum() == 0:
            continue
        ax.plot(cls_data[i + 1]["l"][m], ratio[m],
                color=colors[i], lw=1.2, label=rf"$\beta={b}$")
    ax.axhline(0, color="0.5", lw=0.5)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\ell$")
    ax.set_ylabel(r"$\Delta C_\ell^{TT}/C_\ell^{TT,\Lambda CDM}$")
    ax.set_title("TT")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2, loc="best")

    # --- top-right: Delta C^EE / C^EE ---
    ax = axes[0, 1]
    mask_ee = np.abs(lcdm_cl["EE"]) > 0.05 * np.abs(lcdm_cl["EE"]).max()
    for i, b in enumerate(BETAS):
        ratio = (cls_data[i + 1]["EE"] - lcdm_cl["EE"]) / lcdm_cl["EE"]
        m = mask_ee & finite_only(ratio)
        if m.sum() == 0:
            continue
        ax.plot(cls_data[i + 1]["l"][m], ratio[m],
                color=colors[i], lw=1.2, label=rf"$\beta={b}$")
    ax.axhline(0, color="0.5", lw=0.5)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\ell$")
    ax.set_ylabel(r"$\Delta C_\ell^{EE}/C_\ell^{EE,\Lambda CDM}$")
    ax.set_title("EE")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2, loc="best")

    # --- bottom-left: P(k) absolute curves ---
    ax = axes[1, 0]
    ax.loglog(lcdm_pk["k"], lcdm_pk["P"], color="k", lw=1.5, label="LCDM")
    for i, b in enumerate(BETAS):
        m = finite_only(pks[i + 1]["k"], pks[i + 1]["P"]) & (pks[i + 1]["P"] > 0)
        if m.sum() == 0:
            continue
        ax.loglog(pks[i + 1]["k"][m], pks[i + 1]["P"][m], color=colors[i], lw=1.0,
                  label=rf"$\beta={b}$")
    ax.set_xlabel(r"$k\ [h/\mathrm{Mpc}]$")
    ax.set_ylabel(r"$P(k)\ [(\mathrm{Mpc}/h)^3]$")
    ax.set_title("Linear matter power spectrum")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7, ncol=2)

    # --- bottom-right: Delta P(k) / P(k)^LCDM ---
    ax = axes[1, 1]
    for i, b in enumerate(BETAS):
        # Pks may have slightly different k grids; interpolate to LCDM grid (skip nans)
        good = finite_only(pks[i + 1]["k"], pks[i + 1]["P"]) & (pks[i + 1]["P"] > 0)
        if good.sum() < 2:
            continue
        P_interp = np.interp(lcdm_pk["k"], pks[i + 1]["k"][good], pks[i + 1]["P"][good])
        ratio = (P_interp - lcdm_pk["P"]) / lcdm_pk["P"]
        m = finite_only(ratio)
        if m.sum() == 0:
            continue
        ax.semilogx(lcdm_pk["k"][m], ratio[m], color=colors[i], lw=1.2, label=rf"$\beta={b}$")
    ax.axhline(0, color="0.5", lw=0.5)
    ax.set_xlabel(r"$k\ [h/\mathrm{Mpc}]$")
    ax.set_ylabel(r"$\Delta P(k)/P^{\Lambda CDM}(k)$")
    ax.set_title(r"$\Delta P(k)$ vs $\Lambda CDM$")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7, ncol=2)

    # Set explicit limits to prevent log-tick formatter from hitting log(0)=-inf
    # (matplotlib's auto-layout chokes on it during tight-bbox computation).
    axes[0, 0].set_xlim(2, 2500)
    axes[0, 1].set_xlim(2, 2500)
    axes[1, 0].set_xlim(lcdm_pk["k"].min(), lcdm_pk["k"].max())
    axes[1, 0].set_ylim(lcdm_pk["P"].min() * 0.5, lcdm_pk["P"].max() * 2)
    axes[1, 1].set_xlim(lcdm_pk["k"].min(), lcdm_pk["k"].max())

    fig.suptitle(rf"Forward-model $\beta$ sweep at $c_D={C_D}$, $c_{{s,X}}^{{2}}={CS2_X}$  "
                 r"(Phase 3 pre-flight §5.1)", fontsize=11)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.92, bottom=0.07,
                        wspace=0.25, hspace=0.30)
    fig.savefig(OUTPUT_PDF)
    print(f"\nWrote {OUTPUT_PDF}")

    print("\nNumerical summary (max |Delta C/C| over l masked, plus DeltaP/P at k=0.1 h/Mpc):")
    k_target = 0.1
    for i, b in enumerate(BETAS):
        TT = cls_data[i + 1]["TT"]
        EE = cls_data[i + 1]["EE"]
        if not (np.isfinite(TT).all() and np.isfinite(EE).all() and np.isfinite(pks[i + 1]["P"]).all()):
            print(f"  beta={b:>4}  NaN in CLASS output — pre-flight failure for 3c at this c_D")
            continue
        rTT = np.abs((TT[mask_tt] - lcdm_cl["TT"][mask_tt]) / lcdm_cl["TT"][mask_tt]).max()
        rEE = np.abs((EE[mask_ee] - lcdm_cl["EE"][mask_ee]) / lcdm_cl["EE"][mask_ee]).max()
        P_at_k = np.interp(k_target, pks[i + 1]["k"], pks[i + 1]["P"])
        P_lcdm_at_k = np.interp(k_target, lcdm_pk["k"], lcdm_pk["P"])
        dPP = (P_at_k - P_lcdm_at_k) / P_lcdm_at_k
        print(f"  beta={b:>4}  TT={rTT:.3e}  EE={rEE:.3e}  DeltaP/P|0.1={dPP:+.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
