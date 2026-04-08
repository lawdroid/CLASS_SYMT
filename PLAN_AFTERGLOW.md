# CLASS_SYMT Implementation Plan

**How the Martin & Koh (April 2026) afterglow dark-energy model is layered
onto the public CLASS (Cosmic Linear Anisotropy Solving System) code base.**

This document is the authoritative roadmap for the `CLASS_SYMT` repository.
It explains the *overall approach*, where the new physics plugs into the
existing CLASS module pipeline, which files are touched, and in what order.
It is intentionally written to be readable before any C has been written —
a theorist should be able to audit the plan against the paper without
opening the source tree.

---

## 1. Guiding principles

1. **Additive, not invasive.** CLASS is a tightly-coupled Boltzmann code.
   We add a *new self-contained module* (`afterglow/`) and touch the stock
   modules only at well-defined hook points. When `afterglow_on = 0` the
   run is bit-identical to vanilla CLASS LCDM. This is the same discipline
   used by CLASS's own `dark_radiation`, `fld`, and `scf` add-ons.

2. **Leave the early universe alone.** The Mar 2026 MCMC on the
   `feature/kappa-evolution` branch ruled out any modification of the
   sound horizon (changing Newton's constant shifted the acoustic peaks and
   the other six parameters could not absorb the damage, so the data pulled
   `κ_c` back to 0.998 and `H_0` back to 66.8). The afterglow model
   deliberately deposits all new physics *after* confinement at
   `z_conf ≲ 10⁴`–`10⁶`, leaving `r_s` untouched.

3. **Do not collapse the memory.** The core physical content of the paper
   lives in the Müller–Israel–Stewart memory variable `Σ(τ, x)`. The stock
   CLASS `fld` module algebraically eliminates the analogue of Σ in favour
   of an effective sound speed, which *destroys* the afterglow. We keep
   `δΣ` as an *independent* perturbation DOF at every stage.

4. **Scaffold-first, physics-later.** Phase 0 (this initial commit) is a
   no-op scaffold: the new files compile, `afterglow.ini` runs, and the
   outputs match vanilla CLASS to machine precision. Each subsequent
   phase adds physics on its own branch with a failing-then-passing test.

5. **One branch per phase.** `main` only ever contains code that
   compiles, runs, and matches validation targets. Active development
   happens on `feature/afterglow-*` branches.

---

## 2. Mapping of paper equations to CLASS code

| Paper eq. | Symbol                                      | Where it lives in CLASS_SYMT            |
|-----------|---------------------------------------------|-----------------------------------------|
| Eq. 8, 11 | Yang-Mills Lagrangian, `Λ_D ~ 2.3 meV`      | *external* — sets `z_conf` prior only   |
| Eq. 19    | Imperfect-fluid stress `T^X_μν`             | `source/afterglow/afterglow.c`          |
| Eq. 22    | `p_X = −ρ_X + Σ`                            | `afterglow_background_functions()`      |
| Eq. 30–31 | `ρ_DE = 3 c_D Σ`, `p_DE = −(3 c_D − 1) Σ`   | background index fill                   |
| Eq. 35    | `Ψ(r) = 4r(r−1)/(1+r)³`                     | `afterglow_kernel_Psi()` (implemented)  |
| **Eq. 37**| `u·∇Σ = −(Θ/3c_D)[1 − βΨ(r)] Σ`              | `afterglow_background_derivs()` (stub)  |
| **Eq. 42**| `Q = −(β/c_D) H Ψ(r) ρ_X`                   | `afterglow_exchange_Q()` (implemented)  |
| Eq. 43–44 | Matter-sector backreaction from `Q`         | hook in `background_functions()`        |
| Eq. 46    | `w_eff(r) = −1 + [1 − βΨ(r)] / (3 c_D)`     | derived output only                     |
| Eq. 57    | Linearised continuity+Euler for `X`         | `afterglow_perturb_derivs()`            |
| Eq. 65    | `δΣ` evolution with matter loading source   | `afterglow_perturb_derivs()`            |

Bold rows are the two equations that define the mechanism; everything else
follows by construction once those are in place.

---

## 3. Where afterglow plugs into the CLASS pipeline

CLASS runs its modules in a strict order:

```
input  →  background  →  thermodynamics  →  perturbations  →
          primordial  →  transfer  →  harmonic  →  lensing  →  output
```

The afterglow module hooks into three of them. Nothing downstream of
`perturbations` needs to know the afterglow model exists — it just consumes
a modified `δ_tot(k, τ)` and `H(τ)` like any other cosmology.

### 3.1 `input` (source/input.c)

- Parse the new keys from the `.ini` file (`afterglow_on`, `Omega_dr_dark`,
  `z_conf`, `c_D`, `beta_aft`, `Sigma_today`, `cs2_X`).
- Validate: `c_D > 0`, `beta_aft ≥ 0`, `z_conf > z_rec`, `0 ≤ Omega_dr_dark ≤ 0.1`.
- Copy into `struct afterglow_params` inside `struct precision`.

### 3.2 `background` (source/background.c + source/afterglow/afterglow.c)

- `background_indices()` — reserve indices for `rho_X`, `p_X`, `Sigma`,
  `rho_dr_dark` in the state vector `pvecback`.
- `background_functions()` — called at every `z`; query the afterglow
  module for `ρ_X(z)`, `p_X(z) = −ρ_X + Σ`, and the dark-radiation bath.
  Energy totals (`rho_tot`, `p_tot`) then include them automatically, so
  `H(z)` and the Friedmann solver pick up the new species for free.
- `background_solve_tau()` — add three new lines to the ODE vector:
  `y[rho_X]`, `y[Sigma]`, `y[rho_dr_dark]`, integrated by the existing
  `ndf15` (stiff) or `rkck` (non-stiff) solver.
- Backreaction on CDM (Eq. 43–44): the exchange `Q` is added as a source
  term in the cold-dark-matter continuity equation. Because CLASS evolves
  `ρ_c · a³` analytically (pure `a⁻³`), we switch CDM onto the ODE path
  when `afterglow_on = 1`. This is the one invasive change; guarded by
  the flag so vanilla runs are untouched.

### 3.3 `perturbations` (source/perturbations.c + afterglow module)

- `perturbations_indices()` — register three new perturbation DOFs:
  `delta_X`, `theta_X`, `delta_Sigma`. The last one is the critical
  independent variable — **never algebraically eliminate it**.
- `perturb_derivs()` — at each `k` and `τ`, call
  `afterglow_perturb_derivs()` to fill the RHS for those three entries.
- Gauge: implement in both synchronous and Newtonian. Cross-check the
  ISW tail agrees between the two to < 0.5 %.
- Initial conditions: `δ_X` adiabatic from photon perturbations;
  `δΣ = 0` at `z_conf` (memory seeded only by matter loading).
- Sound speed: rest-frame `cs²_X = 1` by default (causal). Expose as a
  parameter so the prior can later be relaxed if required.

### 3.4 Modules after `perturbations` — no changes

`transfer`, `harmonic`, `lensing`, `output` consume the modified matter
transfer function and `H(z)` with no knowledge of where they came from.

---

## 4. File-by-file change list

Items marked **NEW** are added by CLASS_SYMT. Items marked *(modified)*
are existing CLASS files touched at hook points only.

```
include/afterglow/afterglow.h           NEW    struct + API          [Phase 0 ✔]
source/afterglow/afterglow.c            NEW    kernel, Q, ODE stubs  [Phase 0 ✔]
afterglow.ini                           NEW    sample input          [Phase 0 ✔]
README_AFTERGLOW.md                     NEW    project overview      [Phase 0 ✔]
ATTRIBUTION.md                          NEW    CLASS upstream credit [Phase 0 ✔]
PLAN_AFTERGLOW.md                       NEW    this document         [Phase 0 ✔]

source/input.c                         (mod)   parse new .ini keys   [Phase 1]
include/background.h                   (mod)   indices for X, Σ, dr  [Phase 1]
source/background.c                    (mod)   ODE hooks, bg fns     [Phase 1]
Makefile                               (mod)   compile afterglow/    [Phase 1]

include/perturbations.h                (mod)   3 new pt indices      [Phase 2]
source/perturbations.c                 (mod)   call afterglow RHS    [Phase 2]
source/afterglow/afterglow.c            (ext)   fill perturb_derivs   [Phase 2]

notebooks/CLASS_SYMT_explorer.ipynb     NEW    Colab-ready demo       [Phase 2]
python/mcmc/afterglow_mcmc.py           NEW    emcee driver           [Phase 3]
```

---

## 5. Phased rollout

### Phase 0 — Scaffold (this commit, on `main`)

- Header + stub `.c` compile cleanly.
- `afterglow.ini` with `afterglow_on = 0` runs and produces bit-identical
  output to `default.ini`.
- Documentation (`README_AFTERGLOW.md`, `PLAN_AFTERGLOW.md`,
  `ATTRIBUTION.md`) in place.
- No changes to `background.c` yet.

### Phase 1 — Background sector (~1 week, `feature/afterglow-background`)

**Status (2026-04-09): Phase 1a complete — self-contained integrator passes
all analytic limit tests. Phase 1b (wire into `background.c`) still pending.**

1. **[DONE]** Implement `afterglow_background_rhs()` and
   `afterglow_background_evolve()` as a self-contained 4th-order
   Runge-Kutta integrator in `source/afterglow/afterglow.c`. The three
   DOFs `{ρ_X, Σ, ρ_dr,D}` evolve in `N = ln a` per Eqs. 22, 37, 42.
2. **[DONE]** `test/test_afterglow_bg.py` — pure-Python replica of the
   C RHS + RK4 stepper, verifying analytic limits *without* having to
   compile CLASS. On 2026-04-09 all 17 assertions pass:
   - `β = 0`: `ρ_X(N) = Ω₀_X · e^{-N/c_D}` to 3 × 10⁻⁸ (two values of `c_D`)
   - `c_D → ∞, β = 0`: `ρ_X`, `Σ` constant over 6 e-folds to 6 × 10⁻⁶
   - `|Ψ|_max = 2√3/9` at `r = 2 ± √3`
   - `sign(Ψ)` convention, `Ψ(1) = 0`
   - `sign(Q) = −sign(Ψ)` for `β > 0`
   - numerical smoothness with `β = 2.5, c_D = 1` over 2000 steps
3. **[pending]** Wire `Makefile` + `source/input.c` to parse new `.ini` keys.
4. **[pending]** Add background indices in `background_indices()` and call
   `afterglow_background_evolve()` from `background_solve_tau()`.
5. **[pending]** Feed `Q` into the CDM continuity equation behind the flag.
6. **[pending]** Run full CLASS with `afterglow_on=1` and check
   `r_s` unchanged at < 0.01 % vs LCDM with `z_conf = 10⁵`.
7. Merge to `main` only when the full CLASS integration also passes.

### Phase 2 — Perturbations (~2 weeks, `feature/afterglow-perturbations`)

1. Register `δ_X`, `θ_X`, `δΣ` in `perturbations_indices()`.
2. Implement linearised Eq. 57 + Eq. 65 in synchronous gauge first.
3. Port to Newtonian gauge; cross-check ISW tail (< 0.5 % agreement).
4. Rest-frame sound speed `cs²_X = 1`; expose as parameter.
5. **Validation targets**:
   - TT residual vs LCDM bounded: `max|ΔCℓ/Cℓ| < 5 %` for default
     `(β = 2.5, c_D = 1, z_conf = 10⁵)`.
   - Reduced-Hubble-tension run (`β ≃ 3.0`, `c_D ≃ 1.2`) predicts
     `w_eff(z = 0.5) ≈ −1.03` matching DESI DR2 central value within 1 σ.
   - Disable `δΣ` (set to zero by hand) → model visibly fails the test
     → confirms the memory is doing the work.
6. Ship Colab notebook `CLASS_SYMT_explorer.ipynb` with sliders for
   `β`, `c_D`, `z_conf` and live `Cℓ` + `w_eff(z)` plots.

### Phase 3 — MCMC (~1–2 weeks, `feature/afterglow-mcmc`)

1. Reuse the emcee harness proven on the κ branch (16 walkers, 32 cores,
   Ubuntu workstation).
2. 10-parameter space:
   `{h, ω_b, ω_cdm, n_s, A_s, τ_reio}` + `{z_conf, c_D, β, ΔN_eff}`.
3. Likelihoods: Planck 2018 TTTEEE + lowE + lensing, DESI DR2 BAO,
   Pantheon+ SNe Ia. SH0ES prior in a *separate* tension-aware run.
4. **Success criteria to beat LCDM**:
   - `Δχ² < −9` on Planck + DESI combined (≈ 3 σ preference).
   - `β` posterior excludes zero at ≥ 95 %.
   - `w_eff(z)` posterior tracks DESI DR2 central values within 1 σ.
5. If `c_D → ∞` is posteriorly preferred, report honestly — that means
   the data disfavours the afterglow mechanism at that confidence, just
   as happened with the κ model.

### Phase 4 — Paper integration (~1 week)

1. Generate paper-ready figures directly from `notebooks/`.
2. Write Section 7 (numerical results) of the Martin & Koh paper from
   the MCMC outputs.
3. Tag `v1.0` on `main`.

---

## 6. Data flow at runtime

```
                    ┌──────────────┐
 afterglow.ini  ──▶ │    input.c   │  parse afterglow_params
                    └──────┬───────┘
                           │
                    ┌──────▼───────┐    calls ────▶ afterglow_background_derivs
                    │ background.c │◀── ODE RHS       (Eq. 37, 42)
                    └──────┬───────┘    for X, Σ, dr
                           │   fills rho_tot, p_tot, H(z)
                    ┌──────▼───────┐
                    │thermodynamics│   untouched — uses modified H(z)
                    └──────┬───────┘
                           │
                    ┌──────▼────────┐   calls ────▶ afterglow_perturb_derivs
                    │perturbations.c│◀── 3 new DOFs    (Eq. 57, 65)
                    └──────┬────────┘   δ_X, θ_X, δΣ
                           │
                    ┌──────▼──────┐
                    │ transfer.c  │   untouched
                    └──────┬──────┘
                           │
                    ┌──────▼──────┐
                    │ harmonic.c  │   untouched — produces Cℓ
                    └─────────────┘
```

The diagram makes the "additive, not invasive" principle visible at a
glance: only `input`, `background`, and `perturbations` need to know
about the afterglow, and inside each of those modules the touchpoints
are a handful of labelled hooks.

---

## 7. Risks and mitigations

| Risk                                                  | Mitigation                                                           |
|-------------------------------------------------------|----------------------------------------------------------------------|
| Stiff ODE near `z_conf`                               | Use `ndf15` solver; smooth `Q_c` with `tanh`, width `Δz/z ≈ 0.1`     |
| Gauge artefacts in ISW tail                           | Implement both gauges; freeze `δΣ` super-horizon until `k/aH > 0.1`  |
| `δΣ` accidentally algebraically eliminated            | Code review checklist; unit test that disables `δΣ` and must fail    |
| MCMC finds `c_D → ∞` (LCDM limit preferred)            | Report honestly; this is a valid scientific outcome                  |
| `ΔN_eff` fights tight Planck bound                    | Prior `ΔN_eff ∈ [0, 0.3]`; marginalise `z_conf` if degenerate        |
| OneDrive cannot host `.git` (observed 2026-04-09)     | Keep active git repo outside OneDrive; mirror working tree for users |

---

## 8. Non-goals

- **No Yang-Mills lattice calculation.** `Λ_D ≈ 2.3 meV` and the
  confinement epoch `z_conf` are inputs to CLASS_SYMT, not outputs.
- **No back-compatibility with `fld`.** Users who want the stock fluid
  module keep vanilla CLASS. CLASS_SYMT is focused on one physics target.
- **No modification of recombination or thermodynamics.** `thermodynamics.c`
  is untouched in every phase.
- **No GPU port.** Serial CLASS suffices for MCMC on a 20-core box.

---

## 9. Milestones and acceptance

| Milestone | Artifact                                                      | Done when                                  |
|-----------|---------------------------------------------------------------|---------------------------------------------|
| Phase 0   | Scaffold commit                                               | `afterglow_on=0` matches LCDM bit-for-bit   |
| Phase 1   | Background commit + 4 unit tests                              | `β=0` and `c_D→∞` analytic limits to 10⁻⁶   |
| Phase 2   | Perturbation commit + Colab notebook                          | `w_eff(z)` tracks DESI DR2 within 1 σ       |
| Phase 3   | MCMC corner plot + posterior tables                           | `Δχ² < −9` vs LCDM or honest null report    |
| Phase 4   | Paper Section 7 draft + `v1.0` tag                            | Ready for referee                           |

---

*Document version: 0.1 (Phase 0 scaffold, April 9, 2026).
Updates land on this file with each phase merge.*
