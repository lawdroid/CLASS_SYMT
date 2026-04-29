# Phase 3 MCMC — CLASS_SYMT Afterglow

Three-branch Bayesian fit for the Martin & Koh (April 2026) afterglow
dark-energy model against current cosmological data, run as three sequential
chains. Strategy committed in `Tom/Phase3_MCMC_Master_Plan.md` and
`Tom/Reply_to_Tom_Phase3_Launch.md`. Execution recipe in
`Tom/Phase3_MCMC_Windows_Execution.md`.

## Three-branch design

| Branch | Sampler | β prior | Hypothesis | Runtime | Yaml |
|---|---|---|---|---|---|
| 3a | Cobaya MH × 8 | β = 0 (fixed) | Arrested medium — no DE→matter exchange | 2-4 days / 16 cores | `cobaya_phase3a_baseline.yaml` |
| 3b | Cobaya MH × 8 | β ∈ (0, 2.4] | Non-crossing exchange — DE leaks energy into matter | 7-14 days / 32 cores | `cobaya_phase3b_noncrossing.yaml` |
| 3c | Cobaya MH × 8 | β ∈ [2.65, 3.40] | Crossing strip — effective phantom mimicry | 7-14 days / 32 cores | `cobaya_phase3c_crossing.yaml` |

All three are pre-committed: no dropping 3b or 3c on the basis of 3a results.

## Locked decisions (Tom + Ingyu, 2026-04-30 meeting)

1. ε_Σ via flatness shooting (`afterglow_shoot_eps_Sigma`).
2. `cs2_X = 1/2` fixed across all three branches. Sensitivity at `cs2_X = 1/3`
   verified pre-meeting; see `figures/cb2_sensitivity.pdf`.
3. Late-branch guard `c_D > (1 + 4β)/3` enforced as a hard prior cutoff.
4. Cℓ localization at β = 0.1 verified before 3b launch; see
   `figures/delta_cl_beta_0p1.pdf`.
5. σ₈(z) and f σ₈(z) at z = 0.1, 0.3, 0.5 emitted as derived parameters in
   all three yamls.
6. Pre-registered numerical thresholds in `PHASE3_PREREGISTRATION.md`.

## Datasets

| Probe | Source | Path |
|---|---|---|
| CMB | Planck PR4 (NPIPE) TTTEEE high-l, low-l TT/EE, lensing | `data/Planck_NPIPE/` |
| BAO | DESI DR2 | `data/DESI_DR2/` |
| SNe | Pantheon+ (no SH0ES anchor) | `data/Pantheon+/` |

These directories must be provisioned with the actual data files **before** the
provenance lock (`PHASE3_PREREGISTRATION.md` SHA-256 block) and chain launch.
The submit scripts pre-flight check for their existence.

## Operational targets

| Item | Value |
|---|---|
| Convergence | `Rminus1_stop = 0.01` (tighter than legacy 0.02), ESS > 1000 per parameter |
| Replication | emcee backend, separate from this yaml — to be added |
| Provenance | git tag `phase3-prereg-v1` + data SHA-256 in `PHASE3_PREREGISTRATION.md` |

## Pre-flight diagnostics (already executed)

| Script | Purpose | Output |
|---|---|---|
| `scripts/cb2_sensitivity.py` | Decision 2 — `cs2_X` ∈ {1/3, 1/2} sensitivity at β=0.1 | `figures/cb2_sensitivity.pdf` |
| `scripts/delta_cl_diagnostic.py` | Decision 4 — β=0.1 background-only vs perturbative decomposition | `figures/delta_cl_beta_0p1.pdf` |
| `scripts/forward_model_sweep.py` | Plan §5.1 — Cℓ/Pk at 7 β values across 3a/3b/3c priors | `figures/forward_model_sweep.pdf` |

**Open issue surfaced by `forward_model_sweep`:** at `c_D = 8.0`, the
perturbation integrator produces NaN in P(k) and Cℓ for β ∈ [2.7, 3.4]
(the 3c crossing strip), even though the late-branch guard
`c_D > (1 + 4β)/3 ≈ 4.87` is satisfied. The numerical instability is
narrower than the analytic guard predicts. Resolution required before
3c launch: either narrow the c_D prior in 3c (raise `log10_c_D` floor
above 0.7), or tighten `tol_perturbations_integration` to 1e-7, or
isolate which sub-equation breaks down. Tracked separately as a Tier
follow-up.

## Launch sequence

```powershell
# From CLASS_SYMT/ on Windows PowerShell, after PHASE3_PREREGISTRATION.md
# is committed and the phase3-prereg-v1 tag is pushed:

.\mcmc\submit_phase3a.ps1
# Wait for convergence (Rminus1 < 0.01). Inspect chains/phase3a_baseline/.

.\mcmc\submit_phase3b.ps1
# Same.

.\mcmc\submit_phase3c.ps1   # blocked until the c_D=8.0 NaN issue is resolved
```

Each submit script pre-flights the CLASS binary, dataset directories, and
(for 3b/3c) the prior chain's convergence sentinel.

## Archived files

The earlier two-stage strategy (Stage 1 background-only → Stage 2 PolyChord
perturbations) was superseded by the three-branch strategy on 2026-04-30.
Files renamed `*_archived_2026-04-30.*` are retained for provenance.
