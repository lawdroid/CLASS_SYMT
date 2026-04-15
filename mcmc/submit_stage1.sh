#!/usr/bin/env bash
# ==============================================================
# Stage-1 background MCMC launcher (Linux / macOS / WSL)
# CLASS_SYMT afterglow — Martin & Koh (April 2026)
# Usage:   ./submit_stage1.sh   [from mcmc/ directory]
#          CHAINS=16 ./submit_stage1.sh    (override default)
# ==============================================================
set -euo pipefail

CHAINS=${CHAINS:-8}
YAML="./cobaya_stage1.yaml"
OUT="./chains/stage1"
mkdir -p "$OUT"

# --- sanity checks ---
if [[ ! -x "../class" ]]; then
    echo "[ERROR] CLASS binary not found at ../class."
    echo "        Apply Patches 7-12, then:  cd ..; make class"
    exit 1
fi
if ! command -v cobaya-run >/dev/null 2>&1; then
    echo "[ERROR] Cobaya not installed."
    echo "        Install:  pip install cobaya getdist mpi4py"
    exit 1
fi

cat <<EOF

================================================================
  CLASS_SYMT Phase 3 — Stage 1 (background-only) MCMC
  Datasets: Pantheon+SH0ES · DESI DR2 BAO · Planck 2018 cmp
  Sampler:  Metropolis-Hastings × $CHAINS chains (MPI)
  Target:   R-1 < 0.02 (Gelman-Rubin auto-stop)
================================================================

EOF

# Launch via MPI (preferred) or fall back to parallel processes
if command -v mpirun >/dev/null 2>&1; then
    echo "[INFO] Using MPI (mpirun -n $CHAINS)"
    mpirun -n "$CHAINS" cobaya-run "$YAML" --output-dir "$OUT"
else
    echo "[INFO] mpirun not available — launching $CHAINS serial chains"
    for i in $(seq 1 "$CHAINS"); do
        cobaya-run "$YAML" --output-dir "$OUT/c$i" \
            > "$OUT/chain_$i.log" 2>&1 &
        echo "  launched chain $i (PID $!)"
        sleep 2
    done
    wait
fi

echo ""
echo "Stage 1 complete. Analyze with:  python analyze_stage1.py"
