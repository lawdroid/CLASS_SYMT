# Phase 3 MCMC — CLASS_SYMT Afterglow

Two-stage Bayesian fit for the Martin & Koh (April 2026) afterglow
dark-energy model against current cosmological data.

## Design

| Stage | Sampler | Datasets | Runtime | Purpose |
|-------|---------|----------|---------|---------|
| 1 | Cobaya MCMC × 8 | Pantheon+ · DESI DR2 BAO · Planck 2018 compressed | 2–4 days / 16 cores | Background-only constraint on `c_D`, `β`, `Ω₀_X` |
| 2 | Cobaya PolyChord | Stage-1 + Planck TTTEEE (plik_lite) + lensing + DESI f σ₈ | 7–14 days / 32 cores | Full perturbation fit + Bayes-factor vs ΛCDM |

## Prerequisites

### 1. Apply Patches 7–12 to live CLASS
Before any MCMC can run, you need a CLASS binary that understands the
afterglow parameters. From the repo root:
```powershell
# Apply Patch 7 (Phase 1b.2) to source/background.c
# Apply Patches 8-12 (Phase 2) to include/perturbations.h + source/perturbations.c
# (see PHASE1B2_PATCHES.md and PHASE2_PATCHES.md for exact diffs)

make clean && make class -j8
./class -h | grep afterglow      # should list afterglow_on, afterglow_c_D, afterglow_beta
```

### 2. Install Cobaya + dependencies
```bash
pip install cobaya getdist mpi4py matplotlib
pip install cobaya[polychord]        # Stage-2 nested sampling
cobaya-install cosmo -p ./cobaya_packages   # downloads Planck likelihoods etc
```

### 3. Download datasets
Place under `./data/`:
- `Pantheon+/` — https://github.com/PantheonPlusSH0ES/DataRelease
- `DESI_DR2/` — https://data.desi.lbl.gov/public/dr2
- `Planck2018/` — via `cobaya-install` above

## Run order

### Stage 1 (background-only)
```bash
# Linux / WSL / macOS
./submit_stage1.sh

# Windows PowerShell
.\submit_stage1.ps1
```
Chains auto-stop when Gelman–Rubin R-1 < 0.02. Typical: 60–80k accepted
samples per chain.

### Inspect Stage 1
```bash
python analyze_stage1.py
```
Prints 1D marginals, correlations, Gelman-Rubin, and writes
`chains/stage1/triangle.pdf`.

Headline result to look for:
- `c_D = 1.3 ± 0.15` would match the paper's preferred value.
- `w_X(0) > -1 at 2σ` would be a dynamical-DE detection.

### Stage 2 (full perturbations)
```bash
./submit_stage2.sh
```
PolyChord outputs posterior + evidence log Z. Bayes factor vs ΛCDM
appears in `chains/stage2/afterglow_full.stats`.

## File map

```
mcmc/
├── README.md                     (this file)
├── cobaya_stage1.yaml            Stage-1 config
├── cobaya_stage2.yaml            Stage-2 config
├── afterglow_theory.py           Custom Cobaya theory wrapper
├── submit_stage1.ps1 / .sh       Launchers
├── submit_stage2.sh              Stage-2 launcher
├── analyze_stage1.py             Post-run triage
└── chains/
    ├── stage1/                   Stage-1 chains + triangle.pdf
    └── stage2/                   Stage-2 chains + evidence
```

## Expected outcome (author estimate)

Based on current Planck 2018 + DESI DR2 + Pantheon+ combined analyses
and the paper's fiducial parameter ranges:

- Stage 1 should pin `c_D` to ~10% precision, `β` to ~30%, and
  favor `c_D ≈ 1.3` over `c_D → ∞` (ΛCDM) at 1-2σ if the DESI DR2
  w₀w_a hint holds up in the independent fit.
- Stage 2 will either confirm (tightened c_D, β) or rule out the
  MIS memory signature via growth-sector constraints (σ₈, f σ₈).

Bayes factor log Z_aft − log Z_ΛCDM > 3 would be strong preference
for afterglow; |Δ log Z| < 1 inconclusive; < −3 ruling afterglow
out against the data combo.

## Troubleshooting

**`AfterglowTheory` raises "CLASS does not have afterglow support"**
→ Patches 7–12 not applied or CLASS wasn't rebuilt. Re-run `make class`.

**Chains stuck at R-1 = 0.3+**
→ Check `beta_aft` prior edge — posterior may be pushing against β = 0
boundary. Remove the `min: 0.0` bound or switch to `log10(beta + 1)`.

**PolyChord OOM on 32 cores**
→ Reduce `nlive` from 500 to 300 in `cobaya_stage2.yaml`.

## Citation

Martin & Koh (April 2026), "Dark Energy as the Thermodynamic Afterglow
of a Hidden Gauge-Theory Transition." (Paper attached in `docs/`.)
