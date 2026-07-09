# Phase 3 chains — technical reference

Cobaya MCMC chains from the Martin–Koh Yang-Mills afterglow framework.
This document describes file layout, column structure, parameter
definitions, the likelihood stack, sampler configuration, convergence
diagnostics, and how to read the chains programmatically.

Source repo: https://github.com/ingyukoh/CLASS_SYMT (branch `phase3-interactive`)

---

## 1. Folder layout

```
chains/
    MANIFEST.sha256                          checksums for every file

    phase3a_baseline/                        live, converged
        afterglow_3a.1.txt    135 MB         walker 1 — MCMC samples
        afterglow_3a.2.txt    136 MB         walker 2
        afterglow_3a.3.txt    136 MB         walker 3
        afterglow_3a.4.txt    136 MB         walker 4
        afterglow_3a.checkpoint              sampler state for resume
        afterglow_3a.covmat                  proposal covariance (learnt)
        afterglow_3a.input.yaml              user-supplied input config
        afterglow_3a.updated.yaml            cobaya as-run config (post defaults)
        afterglow_3a.input.yaml.locked       0-byte lock (cobaya housekeeping)
        afterglow_3a.progress                R-1 + acceptance trajectory

    phase3b_noncrossing/                     live, converged (same file pattern)
        afterglow_3b.[1-4].txt    ~210 MB each

    phase3c_b1narrow/                        live, converged (same file pattern)
        afterglow_3c.[1-4].txt    ~340 MB each

    phase3a_baseline_partial_20260506/       HISTORICAL — crashed run (4/30 → 5/06)
                                             R-1 = 0.501 at lensing-Cl crash.
                                             Same file pattern. DO NOT analyze
                                             as final; use phase3a_baseline.

    phase3a_FULL_BACKUP_20260506_1233/       HISTORICAL — backup of the crashed
                                             run including logs/, yaml/, docs/.
                                             Includes README.md describing
                                             the crash mode.
```

The three `live` folders carry the converged chains used for all
published results. The two backup folders are kept for forensics
(why the original 3a crashed at 4 h before reaching R-1 < 0.01).

## 2. Cobaya `.txt` chain format

Each `afterglow_3{a,b,c}.[1-4].txt` is a plain-text MCMC sample file
written by Cobaya. **Line 1** starts with `#` and lists column names,
whitespace-separated. **Lines 2+** are sample rows, one per accepted
proposal, formatted as scientific-notation floats.

**Universal columns (always first two):**
```
weight        integer-valued, count of consecutive rejections after
              this sample (so MCMC density ∝ weight)
minuslogpost  -log(L · π) at this sample, unnormalised
              (i.e. doesn't include log Z)
```

**Next: sampled cosmological + theory params**
```
H0             Hubble constant in km/s/Mpc
omega_b        physical baryon density Ω_b h²
omega_cdm      physical CDM density Ω_cdm h²
tau_reio       reionisation optical depth
logA           ln(10¹⁰ A_s) — primordial amplitude (drop=True; A_s derived)
n_s            scalar spectral index
log10_c_D      log₁₀ of the Yang-Mills coupling c_D (drop=True; c_D derived)
beta_aft       β_aft (only 3b and 3c; in 3a it is fixed at 0 and not in chain)
```

**Next: Planck nuisance params (CamSpec NPIPE TTTEEE)**
```
A_planck       Planck overall calibration
amp_143, amp_217, amp_143x217        extragalactic foreground amplitudes
n_143, n_217, n_143x217              foreground tilts
calTE, calEE                         TE / EE polarisation calibration
```
(Other CamSpec nuisance — `cal0`, `cal2`, `amp_100`, `n_100`,
`use_fg_residual_model` — are fixed at trivial values, see
`input.yaml`, and don't appear in the chain.)

**Next: derived params written by the afterglow theory module + class**
```
A_s                          = 1e-10 * exp(logA)
c_D                          = 10 ** log10_c_D
late_branch_guard            = c_D - (1 + 4β)/3  (regulator diagnostic)
w_X                          = -1 + 1/(3·c_D)    (effective DE EOS)
Omega_m                      late-time matter density
sigma8                       linear growth amplitude at z=0
sigma8_at_z0p{1,3,5}         σ₈ at z = 0.1, 0.3, 0.5
fsigma8_at_z0p{1,3,5}        f·σ₈ at the same z (for RSD)
```

**Last: per-likelihood χ² and total**
```
chi2__BAO, chi2__CMB, chi2__SN
minuslogprior, minuslogprior__0
chi2
chi2__planck_NPIPE_highl_CamSpec.TTTEEE
chi2__planck_2018_lowl.TT
chi2__planck_2018_lowl.EE_sroll2
chi2__planck_2018_lensing.native
chi2__bao.desi_dr2
chi2__sn.pantheonplus
```

To get the header without reading the whole file:
```bash
head -1 chains/phase3a_baseline/afterglow_3a.1.txt
```

## 3. Per-branch sampled / derived / fixed parameter sets

The branches differ ONLY in the treatment of `beta_aft`. Everything
else (priors, likelihoods, nuisance) is identical.

| Param        | 3a baseline | 3b non-crossing | 3c crossing (B1) |
|---|---|---|---|
| `H0`         | sampled U(55, 85)            | sampled U(55, 85)         | sampled U(55, 85)         |
| `omega_b`    | sampled U(0.019, 0.025)      | (same)                    | (same)                    |
| `omega_cdm`  | sampled U(0.1, 0.15)         | (same)                    | (same)                    |
| `tau_reio`   | sampled U(0.01, 0.1)         | (same)                    | (same)                    |
| `logA`       | sampled U(2.5, 3.5)          | (same)                    | (same)                    |
| `n_s`        | sampled U(0.92, 1.0)         | (same)                    | (same)                    |
| `log10_c_D`  | sampled U(0.7, 2.0)          | sampled U(0.7, 2.0)       | sampled U(0.7, 2.0)       |
| `beta_aft`   | **fixed at 0**               | sampled U(0, 2.4)         | sampled U(2.65, 2.75)     |
| Planck nuis  | sampled (see above)          | (same)                    | (same)                    |

Dimensionality of the *sampled* space:
- 3a: d = 7 cosmological + 8 Planck nuisance = **15**
- 3b, 3c: d = 8 cosmological + 8 Planck nuisance = **16**

For the k-NN evidence calculation in `mcmc/evidence_knn.py`, we use
only the 7 (3a) / 8 (3b, 3c) cosmological dimensions — Planck nuisance
marginalizes identically across branches and cancels in Δ log Z.

## 4. Likelihood stack

Identical across all three branches. Specified in each
`afterglow_3*.input.yaml`:

| Likelihood | Type | Source | Notes |
|---|---|---|---|
| `planck_NPIPE_highl_CamSpec.TTTEEE` | CMB | Planck NPIPE CamSpec_NPIPE_12_6 | Uses 143, 217, 143×217, TE, EE |
| `planck_2018_lowl.TT`              | CMB | Planck 2018 PR3 low-ℓ TT (ℓ ∈ [2, 29]) |  |
| `planck_2018_lowl.EE_sroll2`       | CMB | Planck 2018 sroll2 low-ℓ EE |  |
| `planck_2018_lensing.native`       | CMB lensing | Planck 2018 SMICA DX12 lensing | mv2 ndclpp_p teb_consext8 |
| `bao.desi_dr2`                     | BAO | DESI DR2 — GCcomb mean + cov | |
| `sn.pantheonplus`                  | SN Ia | Pantheon+ — use_abs_mag: false | M_B marginalized |

All four CMB likelihoods carry `stop_at_error: false` (cobaya-block
level + per-likelihood level — see lesson §3.12 in CLAUDE.md) so a
single missing-Cl exception doesn't abort the chain.

## 5. Theory module

`afterglow_theory.AfterglowTheory`, a thin wrapper around CLASS that
adds the YM afterglow extension. Key extra_args:

```yaml
output: tCl,pCl,lCl,mPk
l_max_scalars: 2508
lensing: 'yes'
non_linear: halofit
afterglow_on: 1            # turn the YM extension on
cs2_X: 0.5                 # dark-medium sound speed (fixed)
tol_background_integration: 1.0e-05
tol_perturbations_integration: 1.0e-07
z_max_pk: 1.0
```

Version `v3.3.4` (`afterglow_3*.updated.yaml`). The `late_branch_guard`
output column is the smooth regulator added 2026-05-15 to avoid
sigma_hat sign-flip near c_D = (1 + 4β)/3 — see CLAUDE.md §3.12 and
`PHASE1B_PATCHES.md` in the repo.

## 6. Sampler configuration

Cobaya's adaptive Metropolis-Hastings MCMC, four parallel walkers
across MPI ranks. From `updated.yaml`:

```yaml
sampler:
  mcmc:
    Rminus1_stop: 0.01              # convergence target
    Rminus1_cl_stop: 0.05
    proposal_scale: 2.4
    learn_proposal: true
    learn_proposal_Rminus1_max: 2.0
    learn_proposal_Rminus1_max_early: 30.0
    oversample_power: 0.4           # fast/slow blocking
    oversample_thin: true
    blocking:
      - [1, [H0, omega_b, omega_cdm, tau_reio, logA, n_s, log10_c_D]]
      - [7, [A_planck]]
      - [7, [amp_143, amp_217, amp_143x217, n_143, n_217, n_143x217,
             calTE, calEE]]
    max_samples: 500000
    burn_in: 0                      # discarded post-hoc in analysis
```

The `blocking` defines fast/slow parameter groups: cosmological params
recompute the full Boltzmann hierarchy (slow); Planck nuisance only
re-evaluates the likelihood (fast, oversampled 7×).

## 7. Convergence diagnostics

Gelman-Rubin R−1 trajectory is in `afterglow_3*.progress` — a TSV with
columns `N, timestamp, acceptance_rate, Rminus1, Rminus1_cl`. Tail of
the 3a progress file at end of run:

```
N         time                 acceptance   R−1        R−1_cl
715949    2026-06-02 06:28:50  0.7040       0.0078     0.0552
718692    2026-06-02 08:45:14  0.7039       0.0076     0.0539
721381    2026-06-02 10:52:19  0.7039       0.0077     0.0547
```

Final values per branch (from `mcmc/export_chains_summary.py`):

| Branch | n_samples (weighted) | n_eff | R−1 |
|---|---|---|---|
| 3a | 505,327 | 387k | 0.0013 |
| 3b | 764,257 | 581k | 0.0002 |
| 3c | 1,250,840 | 1.01M | 0.0004 |

All chains pass R−1 < 0.01 by a factor of ~100×, which is well past
the field-standard "converged" threshold.

## 8. Reading the chains in Python

### Quick numpy (loads all samples for one walker)

```python
import numpy as np
data = np.genfromtxt('chains/phase3a_baseline/afterglow_3a.1.txt',
                      names=True, deletechars='', skip_header=0)
# data['weight'], data['minuslogpost'], data['H0'], data['c_D'], ...
```

### Concatenate all four walkers with burn-in cut

```python
import numpy as np, glob, os
def load_branch(branch_dir, burn_frac=0.3):
    files = sorted(glob.glob(os.path.join(branch_dir, '*.[1-9].txt')))
    parts = []
    for f in files:
        d = np.genfromtxt(f, names=True, deletechars='')
        cut = int(len(d) * burn_frac)
        parts.append(d[cut:])
    return np.concatenate(parts)
samples = load_branch('chains/phase3a_baseline')
```

### Recommended: getdist (handles weights, plots posteriors)

```python
from getdist import loadMCSamples
samples = loadMCSamples('chains/phase3a_baseline/afterglow_3a',
                        settings={'ignore_rows': 0.3})
print(samples.getInlineLatex('c_D', limit=2))   # 95% c_D bound
samples.getCorrelationMatrix()
```

### Effective sample size and per-param marginals

```python
w = samples.weights
n_eff = w.sum()**2 / (w**2).sum()
print(f'n_eff = {n_eff:.0f}')
# Marginal on beta_aft (3b chain):
samples = loadMCSamples('chains/phase3b_noncrossing/afterglow_3b',
                        settings={'ignore_rows': 0.3})
print(samples.getInlineLatex('beta_aft'))    # will report "flat" / wide bounds
```

### Recomputing log Z (Bayesian evidence)

```bash
# In the repo root:
python mcmc/evidence_knn.py
# Reads chains/phase3{a,b,c}_*/ and writes evidence_knn_results.json
# Method: Heavens 2017 k-NN density estimator (arXiv:1704.03472)
```

## 9. Historical backup folders

`phase3a_baseline_partial_20260506/` and
`phase3a_FULL_BACKUP_20260506_1233/` are kept solely for forensics.
The original 3a run crashed at 9 h 45 min wall time with R−1 = 0.501
when `planck_2018_lensing.native` raised "No lensed Cl's were
computed" and aborted MPI. The fix (Phase 3 / B1
smooth regulator + per-likelihood `stop_at_error: false`) is in the
live `phase3a_baseline/` chain.

If you want to inspect the failure: the FULL_BACKUP includes
`logs/<latest>.log` with the abort trace, `chain/` with the partial
samples at R−1 = 0.5, and `docs/` with my hand-written diagnosis.

## 10. Reproducibility

```bash
git clone https://github.com/ingyukoh/CLASS_SYMT.git
cd CLASS_SYMT
git checkout phase3-interactive

# To verify the on-disk chains against MANIFEST.sha256:
cd chains
sha256sum -c MANIFEST.sha256        # Linux/Mac
# or on Windows:
Get-FileHash * -Algorithm SHA256

# To re-launch any branch (Ubuntu workstation pattern):
cobaya-run mcmc/cobaya_phase3a_baseline.yaml --resume
cobaya-run mcmc/cobaya_phase3b_noncrossing.yaml --resume
cobaya-run mcmc/cobaya_phase3c_b1narrow.yaml --resume

# To regenerate the chains_summary.json + k-NN evidence:
python mcmc/export_chains_summary.py
python mcmc/evidence_knn.py
```

Cobaya version `3.6.2`; CLASS version derived from the `theory:` block
of the input yaml. Python ≥ 3.9 expected.

## 11. Contact and provenance

Data: Ingyu Koh
Repo: https://github.com/lawdroid/CLASS_SYMT
Generated 2026-06-02 / 2026-06-03 from chains running on:
- Ubuntu workstation (3a, 3b)
- Mac mini M4 (3c)
