#!/usr/bin/env python3
"""
cb2_sensitivity.py — c_b^2 (cs2_X) sensitivity diagnostic for Phase 3 decision 2.

Runs CLASS_SYMT three times at fixed afterglow background (c_D = 8.0, beta = 0.1):
    - cs2_X = 1/2 (Tom's proposed default for 3a)
    - cs2_X = 1/3 (Ingyu's sensitivity probe per Reply email)
    - LCDM reference (afterglow_on = 0)

Plots Delta C_l / C_l vs l for TT / EE / TE on a 1x3 grid, saved to
figures/cb2_sensitivity.pdf. The visible separation between the two cs2_X
curves answers "is the c_b^2 = 1/2 choice amplifying the beta = 0.1 shift,
or is the result closure-independent?" — the question Ingyu raised in
Reply_to_Tom_Phase3_Launch.md decision 2.

Mapping caveat: the CLASS source calls this parameter `cs2_X`. The Reply
email refers to it as `c_b^2`. Whether these are the same quantity in the
paper's notation needs Tom's confirmation in the 4/30 meeting.

Usage (from CLASS_SYMT root):
    python3 scripts/cb2_sensitivity.py
"""
from __future__ import annotations
import os
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
WORK_DIR = Path("/tmp/cb2_sensitivity_work")
FIGURES_DIR = CLASS_ROOT / "figures"
OUTPUT_PDF = FIGURES_DIR / "cb2_sensitivity.pdf"

# Fixed cosmology (Planck 2018 best-fit-ish, doesn't matter for the diagnostic)
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
    {"label": "LCDM",          "afterglow_on": 0, "c_D": None, "beta": None, "cs2_X": None,         "color": "k",       "ls": "-"},
    {"label": r"$c_{s,X}^{2}=1/2$", "afterglow_on": 1, "c_D": 8.0,  "beta": 0.1,  "cs2_X": 0.5,           "color": "C0",      "ls": "-"},
    {"label": r"$c_{s,X}^{2}=1/3$", "afterglow_on": 1, "c_D": 8.0,  "beta": 0.1,  "cs2_X": 1.0 / 3.0,    "color": "C1",      "ls": "--"},
]


# ---------- helpers ----------
def write_ini(case_idx: int, case: dict) -> Path:
    """Generate the .ini file for a single case and return its path."""
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
    """Invoke ./class on the .ini and raise if it fails."""
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
    """Read CLASS lensed Cl output. Returns dict with l, TT, EE, TE arrays."""
    path = WORK_DIR / f"case_{case_idx:02d}_00_cl_lensed.dat"
    arr = np.loadtxt(path)
    return {"l": arr[:, 0], "TT": arr[:, 1], "EE": arr[:, 2], "TE": arr[:, 3]}


# ---------- main ----------
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
    titles = [r"$\Delta C_\ell^{TT}/C_\ell^{TT}$",
              r"$\Delta C_\ell^{EE}/C_\ell^{EE}$",
              r"$\Delta C_\ell^{TE}/C_\ell^{TE}$"]

    for ax, spec, title in zip(axes, spectra, titles):
        # Mask zero-crossings (mostly relevant for TE which crosses zero, blowing up the ratio)
        mask = np.abs(lcdm[spec]) > 0.05 * np.abs(lcdm[spec]).max()
        for i in (1, 2):  # skip LCDM itself
            case = CASES[i]
            ratio = (cls_data[i][spec] - lcdm[spec]) / lcdm[spec]
            ax.plot(cls_data[i]["l"][mask], ratio[mask],
                    label=case["label"], color=case["color"], ls=case["ls"], lw=1.4)
        ax.axhline(0, color="0.5", lw=0.5)
        ax.set_xlabel(r"$\ell$")
        ax.set_xscale("log")
        ax.set_title(title)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)

    axes[0].set_ylabel(r"$\Delta C_\ell / C_\ell^{\Lambda CDM}$")
    fig.suptitle(r"$c_{s,X}^{2}$ sensitivity at $\beta=0.1$, $c_D=8.0$  (Phase 3 decision 2 probe)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(OUTPUT_PDF, bbox_inches="tight")
    print(f"\nWrote {OUTPUT_PDF}")

    # Print a one-line numerical summary so the result is visible without opening the PDF.
    # Apply the same zero-crossing mask used for plotting (relevant for TE).
    print("\nNumerical summary (max |Delta C/C| over ell, masking |C^LCDM| < 5% peak):")
    for spec in spectra:
        mask = np.abs(lcdm[spec]) > 0.05 * np.abs(lcdm[spec]).max()
        r1 = np.abs((cls_data[1][spec][mask] - lcdm[spec][mask]) / lcdm[spec][mask])
        r2 = np.abs((cls_data[2][spec][mask] - lcdm[spec][mask]) / lcdm[spec][mask])
        print(f"  {spec}:  cs2=1/2  max={r1.max():.3e}   cs2=1/3  max={r2.max():.3e}   spread={abs(r1.max()-r2.max()):.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
