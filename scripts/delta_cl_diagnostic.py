#!/usr/bin/env python3
"""
delta_cl_diagnostic.py — Cl localization diagnostic for Phase 3 decision 4.

Decomposes the beta = 0.1 afterglow effect on the CMB into background-only
vs perturbative contributions. Three CLASS_SYMT runs at fixed c_D = 8.0,
cs2_X = 0.5:
    - LCDM reference (afterglow_on = 0)
    - Afterglow with beta = 0   (background-only difference vs LCDM)
    - Afterglow with beta = 0.1 (full physics)

Plots Delta C_l / C_l^LCDM vs l for TT / EE / TE on a 1x3 grid, with both
afterglow curves overlaid. The vertical separation between the two
afterglow curves at any l is the *perturbative* contribution from the
beta-dependent dark-sector exchange. The shared offset is the background-
only effect from the modified H(z).

Per Reply_to_Tom_Phase3_Launch.md decision 4, this is the pre-3b sanity
check on the "300% Cl shift" magnitude: if the perturbative piece is
>10x off from analytic intuition, it's an implementation bug, not physics.

Usage (from CLASS_SYMT root):
    python3 scripts/delta_cl_diagnostic.py
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
WORK_DIR = Path("/tmp/delta_cl_diagnostic_work")
FIGURES_DIR = CLASS_ROOT / "figures"
OUTPUT_PDF = FIGURES_DIR / "delta_cl_beta_0p1.pdf"

BASE_INI = """
h = 0.6732
omega_b = 0.02238
omega_cdm = 0.1201
A_s = 2.1e-9
n_s = 0.966
tau_reio = 0.0543
output = tCl,pCl,lCl
lensing = yes
l_max_scalars = 2500
"""

CASES = [
    {"label": "LCDM",                 "afterglow_on": 0, "c_D": None, "beta": None, "cs2_X": None, "color": "k",  "ls": "-"},
    {"label": r"AG, $\beta=0$",       "afterglow_on": 1, "c_D": 8.0,  "beta": 0.0,  "cs2_X": 0.5,  "color": "C2", "ls": "--"},
    {"label": r"AG, $\beta=0.1$",     "afterglow_on": 1, "c_D": 8.0,  "beta": 0.1,  "cs2_X": 0.5,  "color": "C3", "ls": "-"},
]


def write_ini(case_idx: int, case: dict) -> Path:
    root = WORK_DIR / f"case_{case_idx:02d}_"
    ini_lines = [BASE_INI, f"root = {root}"]
    if case["afterglow_on"] == 1:
        ini_lines += [
            "afterglow_on = 1",
            f"c_D = {case['c_D']}",
            f"beta_aft = {case['beta']}",
            f"cs2_X = {case['cs2_X']}",
        ]
    ini_path = WORK_DIR / f"case_{case_idx:02d}.ini"
    ini_path.write_text("\n".join(ini_lines) + "\n")
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
        raise RuntimeError(f"CLASS failed on {ini_path} (exit {proc.returncode})")


def read_cl_lensed(case_idx: int) -> dict:
    path = WORK_DIR / f"case_{case_idx:02d}_00_cl_lensed.dat"
    arr = np.loadtxt(path)
    return {"l": arr[:, 0], "TT": arr[:, 1], "EE": arr[:, 2], "TE": arr[:, 3]}


def main() -> int:
    if not CLASS_BIN.exists():
        sys.stderr.write(f"ERROR: CLASS binary not found at {CLASS_BIN}. Run `make class` first.\n")
        return 1

    if WORK_DIR.exists():
        shutil.rmtree(WORK_DIR)
    WORK_DIR.mkdir(parents=True)
    FIGURES_DIR.mkdir(exist_ok=True)

    print(f"Running {len(CASES)} CLASS configurations:")
    for i, case in enumerate(CASES):
        print(f"  case {i}: {case['label']}")
        ini = write_ini(i, case)
        run_class(ini)

    cls_data = [read_cl_lensed(i) for i in range(len(CASES))]
    lcdm = cls_data[0]

    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharex=True)
    spectra = ["TT", "EE", "TE"]
    titles = [r"$\Delta C_\ell^{TT}/C_\ell^{TT,\Lambda CDM}$",
              r"$\Delta C_\ell^{EE}/C_\ell^{EE,\Lambda CDM}$",
              r"$\Delta C_\ell^{TE}/C_\ell^{TE,\Lambda CDM}$"]

    for ax, spec, title in zip(axes, spectra, titles):
        mask = np.abs(lcdm[spec]) > 0.05 * np.abs(lcdm[spec]).max()
        for i in (1, 2):
            case = CASES[i]
            ratio = (cls_data[i][spec] - lcdm[spec]) / lcdm[spec]
            ax.plot(cls_data[i]["l"][mask], ratio[mask],
                    label=case["label"], color=case["color"], ls=case["ls"], lw=1.4)
        # Also show the perturbative-only contribution: (beta=0.1 - beta=0)
        pert = (cls_data[2][spec] - cls_data[1][spec]) / lcdm[spec]
        ax.plot(cls_data[2]["l"][mask], pert[mask],
                label=r"perturbative ($\Delta\beta$)", color="C5", ls=":", lw=1.2)
        ax.axhline(0, color="0.5", lw=0.5)
        ax.set_xlabel(r"$\ell$")
        ax.set_xscale("log")
        ax.set_title(title)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)

    axes[0].set_ylabel(r"$\Delta C_\ell / C_\ell^{\Lambda CDM}$")
    fig.suptitle(r"Afterglow vs $\Lambda CDM$ at $c_D=8.0$, $c_{s,X}^{2}=0.5$  "
                 r"(Phase 3 decision 4: $\beta=0.1$ Cl localization)", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUTPUT_PDF, bbox_inches="tight")
    print(f"\nWrote {OUTPUT_PDF}")

    print("\nNumerical summary (max |Delta C/C^LCDM| over ell, masked):")
    for spec in spectra:
        mask = np.abs(lcdm[spec]) > 0.05 * np.abs(lcdm[spec]).max()
        rb0 = np.abs((cls_data[1][spec][mask] - lcdm[spec][mask]) / lcdm[spec][mask])
        rb1 = np.abs((cls_data[2][spec][mask] - lcdm[spec][mask]) / lcdm[spec][mask])
        rpert = np.abs((cls_data[2][spec][mask] - cls_data[1][spec][mask]) / lcdm[spec][mask])
        print(f"  {spec}:  bg-only(beta=0)={rb0.max():.3e}   full(beta=0.1)={rb1.max():.3e}   pert-only={rpert.max():.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
