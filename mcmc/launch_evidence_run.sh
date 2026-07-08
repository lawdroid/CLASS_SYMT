#!/usr/bin/env bash
# ==============================================================
# Definitive Bayesian-evidence run for  Delta log Z (afterglow - LCDM).
#
# Runs PolyChord nested sampling for BOTH models on the IDENTICAL data
# stack with the SAME sampler, so the evidence difference is exact:
#   1. LCDM reference          -> chains/evidence_lcdm/
#   2. Afterglow 3a (beta=0)   -> chains/evidence_3a/   (nested model)
#
# The Savage-Dickey pre-estimate (evidence_vs_lcdm_sddr.py) gives
# Delta log Z(3a-LCDM) ~ +0.5; this run replaces the boundary
# extrapolation with a direct, absolute log Z for each model.
#
# Run on the compute workstation:
#   cd <CLASS_SYMT_ROOT>
#   nohup bash mcmc/launch_evidence_run.sh > evidence_run.log 2>&1 &
#
# Expected wall time: ~8-12 hr each on 16 MPI ranks (~20 sampled dims,
# nlive=500, num_repeats~3*nDim).  Run sequentially or on two boxes.
# ==============================================================
set -euo pipefail
cd "$(dirname "$0")/.."          # CLASS_SYMT root

NPROC="${NPROC:-16}"
COBAYA_RUN="${COBAYA_RUN:-cobaya-run}"

echo "[$(date)] Evidence run start; NPROC=$NPROC"

mkdir -p chains/evidence_lcdm chains/evidence_3a

# --- 1. LCDM reference (stock classy) ---
echo "[$(date)] === LCDM PolyChord ==="
mpirun -np "$NPROC" "$COBAYA_RUN" mcmc/cobaya_lcdm_polychord.yaml --resume \
  || mpirun -np "$NPROC" "$COBAYA_RUN" mcmc/cobaya_lcdm_polychord.yaml

# --- 2. Afterglow 3a (beta=0), same sampler ---
echo "[$(date)] === Afterglow-3a PolyChord ==="
mpirun -np "$NPROC" "$COBAYA_RUN" mcmc/cobaya_phase3a_polychord.yaml --resume \
  || mpirun -np "$NPROC" "$COBAYA_RUN" mcmc/cobaya_phase3a_polychord.yaml

# --- 3. Read off the evidences and report Delta log Z ---
echo "[$(date)] === Delta log Z ==="
python3 - <<'PY'
import re, glob, sys
def logZ(stats_glob):
    f = sorted(glob.glob(stats_glob))
    if not f:
        print(f"  MISSING: {stats_glob}"); return None
    txt = open(f[0]).read()
    m = re.search(r"log\(Z\)\s*=\s*([-\d.]+)\s*\+/-\s*([\d.]+)", txt)
    if not m:
        m = re.search(r"Z\s*=.*?log.*?=\s*([-\d.]+).*?([\d.]+)", txt)
    return (float(m.group(1)), float(m.group(2))) if m else None
lc = logZ("chains/evidence_lcdm/*.stats")
af = logZ("chains/evidence_3a/*.stats")
print("LCDM   log Z =", lc)
print("3a     log Z =", af)
if lc and af:
    import math
    d = af[0]-lc[0]; e = math.hypot(af[1], lc[1])
    print(f"Delta log Z (3a - LCDM) = {d:+.3f} +/- {e:.3f}")
    print("(+ favors afterglow; Jeffreys: |d|<1 inconclusive, 1-2.5 weak, 2.5-5 moderate, >5 strong)")
PY
echo "[$(date)] Evidence run done."
