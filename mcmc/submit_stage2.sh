#!/usr/bin/env bash
# ==============================================================
# Stage-2 full perturbations MCMC launcher
# CLASS_SYMT afterglow — Martin & Koh (April 2026)
#
# Uses PolyChord (nested sampling): gives evidence Z in addition
# to posterior, so Stage-2 output can do Bayes-factor comparison
# vs LCDM in one pass.
#
# Usage:   ./submit_stage2.sh   [from mcmc/ directory]
# ==============================================================
set -euo pipefail

YAML="./cobaya_stage2.yaml"
OUT="./chains/stage2"
STAGE1_COVMAT="./chains/stage1/afterglow_bg.covmat"
NCORES=${NCORES:-32}

mkdir -p "$OUT"

# --- prereq checks ---
if [[ ! -f "$STAGE1_COVMAT" ]]; then
    echo "[WARN] Stage-1 covmat not found at $STAGE1_COVMAT"
    echo "       PolyChord will work without it but start less efficient."
fi
if ! python -c "import pypolychord" 2>/dev/null; then
    echo "[ERROR] PolyChord not installed."
    echo "        Install:  pip install cobaya[polychord]"
    exit 1
fi

cat <<EOF

================================================================
  CLASS_SYMT Phase 3 — Stage 2 (full perturbations) MCMC
  Datasets: Stage-1 + Planck plik_lite + lensing + DESI f sigma_8
  Sampler:  PolyChord nested (nlive=500, evidence enabled)
  Cores:    $NCORES (MPI)
================================================================

EOF

mpirun -n "$NCORES" cobaya-run "$YAML" --output-dir "$OUT"

echo ""
echo "Stage 2 complete. Analyze with:  python analyze_stage2.py"
echo "Bayes factor LCDM vs afterglow: see chains/stage2/afterglow_full.stats"
