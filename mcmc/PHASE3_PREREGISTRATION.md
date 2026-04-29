# Phase 3 Pre-Registration

**Date:** 2026-04-30
**Authors:** Ingyu Koh, Tom Martin
**Code commit at registration:** *to be filled in by the
`phase3-prereg-v1` git tag at the moment the data SHA-256 block below is
populated and the chains are launched.*
**Paper:** Martin & Koh (April 2026), "Dark Energy as the Thermodynamic
Afterglow of a Hidden Gauge-Theory Transition."

This document fixes, before any chain output is examined, the numerical
thresholds against which the three Phase 3 branches will be judged. After
the `phase3-prereg-v1` tag is signed and pushed, the predictions below are
not editable. Any change becomes a `phase3-prereg-v2` tag with rationale
documented here.

---

## 1. Datasets and provenance lock

Three data products are used. SHA-256 checksums of every file consumed by
the chains will be filled in below before launch.

| Probe | Source | Path | SHA-256 |
|---|---|---|---|
| CMB | Planck PR4 (NPIPE) TTTEEE high-l + low-l TT/EE + lensing | `data/Planck_NPIPE/<files>` | *to be filled* |
| BAO | DESI DR2 | `data/DESI_DR2/<files>` | *to be filled* |
| SNe | Pantheon+ (no SH0ES anchor) | `data/Pantheon+/<files>` | *to be filled* |

The fill-in step is part of plan §7 (provenance lock) and happens after
data files are downloaded but before the first `submit_phase3*.ps1` run.

## 2. Locked decisions (Tom + Ingyu, 2026-04-30 meeting)

These six choices are locked per `Tom/Reply_to_Tom_Phase3_Launch.md`:

1. ε_Σ via flatness shooting (`afterglow_shoot_eps_Sigma`), all three
   branches.
2. `cs2_X = 1/2` fixed across all three branches. Sensitivity at
   `cs2_X = 1/3` verified pre-meeting (`figures/cb2_sensitivity.pdf`):
   TT spread ~0.9%, EE/TE closure-independent.
3. Late-branch guard `c_D > (1 + 4 β) / 3` enforced as a hard prior
   cutoff in 3b and 3c (vacuous in 3a at β=0).
4. β=0.1 Cℓ localization verified pre-meeting
   (`figures/delta_cl_beta_0p1.pdf`): perturbative-only contribution is
   ~0.78% of LCDM in TT, ~0.28% in EE — much smaller than the "300%
   shift" concern would imply.
5. σ₈ and f σ₈ at z = 0.1, 0.3, 0.5 emitted as derived parameters in
   all three yamls.
6. Pre-registered numerical thresholds per §3 below.

## 3. Pre-registered numerical predictions

### Branch 3a (β = 0, baseline / arrested-medium hypothesis)

- **Δχ² ≲ 4** vs ΛCDM on the combined Planck + DESI + Pantheon+
  likelihood. Larger Δχ² would falsify the arrested-medium hypothesis
  at the data-combo level.
- **c_D posterior concentrates at c_D > 5** (i.e. log₁₀ c_D > 0.7).
  A peak at c_D ≈ 1 — the paper's canonical benchmark — would
  require revisiting the prior choice but is itself the *publishable*
  signal.
- **Intrinsic w_X = −1 + 1/(3 c_D) in [−0.98, −0.95]** at the c_D
  posterior mean. w_X < −1 at any sample is forbidden classically and
  would indicate a numerical bug, not new physics.

### Branch 3b (β ∈ (0, 2.4], non-crossing strip)

- **|Δ σ₈| ≥ 0.01** across the β range, in the direction relieving
  the S₈ tension (lower σ₈ relative to ΛCDM at the BAO+CMB-preferred
  Ω_m). A null result at this level documents that the exchange term
  is observationally inert in this strip — itself a publishable
  contribution.
- **Δχ² between 3a (β=0 fixed) and 3b (β free)** to be reported but
  not pre-thresholded; will inform whether 3b adds explanatory power.

### Branch 3c (β ∈ [2.65, 3.40], crossing strip / phantom mimicry)

- **w_eff dips to ≈ −1 − 10⁻²** in z ∈ [0.5, 2.0], computed from the
  posterior mean trajectory. Larger phantom-side excursions indicate
  numerical instability rather than new physics.
- **Intrinsic w_X stays ≥ −1** at every sample (hard floor enforced
  by `c_D > 1/3` part of the late-branch guard). Any sample with
  w_X < −1 indicates an integrator failure.

### Microscopic reconstruction (post-3a/b/c)

The chain over (z_c, ΔN_eff, c_D) maps to the microscopic
parameter ξ via the dimensional-transmutation relation in
`Tom/Phase3_MCMC_Master_Plan.md`. Pre-registered intervals:

- **ξ posterior in [0, 0.30]** (consistent with cosmological
  constraints).
- **ξ < 0.15:** Simons Observatory will not detect the dark gluon
  bath at design sensitivity — null result for SO.
- **ξ ∈ [0.15, 0.30]:** SO detects at ≥ 2σ — confirmation channel.
- **ξ > 0.30:** already excluded by current constraints; would
  indicate inconsistency between the chain and existing bounds.

## 4. Pre-flight findings (already executed)

Three pre-flight scripts were run before the meeting; results are in
`figures/`:

| Script | Result |
|---|---|
| `cb2_sensitivity.py` | TT shift 3.6-4.4% at β=0.1; cs2 spread <0.9%; EE/TE closure-independent. |
| `delta_cl_diagnostic.py` | Background-only contribution dominates (TT 3.4%, EE 5.9%); perturbative piece is ~5% of total. |
| `forward_model_sweep.py` | Smooth Cℓ/Pk evolution for β ∈ [0, 2.0]; **NaN/Inf for β ≥ 2.7 at c_D = 8.0** — see §5. |

## 5. Open issue blocking 3c launch

The forward-model sweep (`figures/forward_model_sweep.pdf`) shows that
at `c_D = 8.0` and `cs2_X = 0.5`, the CLASS_SYMT perturbation
integrator produces non-finite output (NaN in P(k) and Cℓ) for
β ∈ {2.7, 3.0, 3.4}. The late-branch guard `c_D > (1 + 4β)/3`
requires only c_D > 4.87 at β = 3.4 and is satisfied; the practical
numerical-stability boundary is narrower than the analytic guard
predicts.

**Resolution required before 3c launch:**

1. Investigate which sub-equation diverges at β ≥ 2.65 — likely
   candidate is the σ̂ memory variable equation (P3 in
   `source/afterglow/afterglow_pert.c`) hitting a stiff regime when
   β Ψ approaches unity.
2. Tighten `tol_perturbations_integration` from 1e-5 to 1e-7 and
   re-run the sweep at β = 2.7. If that fixes it, lock the tighter
   tolerance into `cobaya_phase3c_crossing.yaml`.
3. If the issue persists, raise the `log10_c_D` floor in the 3c yaml
   from 0.7 (c_D > 5) to ~1.0 (c_D > 10) and document the
   prior-narrowing.

3a and 3b are unaffected by this issue (β = 0 and β ≤ 2.4 are stable
in the sweep) and may proceed on schedule.

## 6. Operational targets

| Target | Value |
|---|---|
| Convergence | `Rminus1_stop = 0.01`, ESS > 1000 per parameter |
| Replication | emcee backend per branch, posterior consistency check |
| Logging | `output_every: 60s`, `chains/<branch>/afterglow_<branch>.{1,2,...}.txt` |
| Re-run policy | Chain extension on `Rminus1 ≥ 0.05` after `max_samples` allowed without re-tag |

## 7. Tag and freeze

Once `data/<probe>/` directories are populated and the SHA-256 block in
§1 is filled, the registration is frozen by:

```
git add mcmc/PHASE3_PREREGISTRATION.md mcmc/cobaya_phase3*.yaml
git commit -m "Phase 3 provenance lock: data SHA-256 + final yamls"
git tag -a phase3-prereg-v1 -m "Phase 3 pre-registration, frozen 2026-XX-XX"
git push origin main
git push origin phase3-prereg-v1
git push ingyukoh main
git push ingyukoh phase3-prereg-v1
```

After the tag exists, no edits to this file or to the three yamls
without a v2 tag and a rationale paragraph appended below.

## 8. Revision history

- **v1 (this file, 2026-04-30):** initial draft. Datasets identified,
  yamls drafted, three branches pre-registered, pre-flight executed.
  Open issue: 3c numerical stability at c_D=8.0 (§5). Awaiting data
  provisioning + SHA-256 fill-in.
