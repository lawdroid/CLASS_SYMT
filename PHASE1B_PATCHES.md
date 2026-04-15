# Phase 1b — Wire the afterglow fluid into CLASS's background.c

Phase 1a delivered a self-contained background integrator validated to
1e-9 against four analytic limits of Martin & Koh (April 2026).
Phase 1b wires that integrator into CLASS's own background solver by
following the `has_fld` pattern at four insertion points in
`background.c`, plus one new block in `input.c` for `.ini` parsing.

All new physics stays inside `source/afterglow/afterglow_class_glue.c`
and `include/afterglow/afterglow_class_glue.h` — `background.c` only
ever calls three glue helpers whose semantics are fully pinned by
Tests 7, 8, 9 in `test/test_afterglow_bg.py`.

---

## Paper equations used (all VERIFIED against `docs/Martin_Koh_Afterglow_April2026.pdf`)

| Eq. | Section | Role in Phase 1b |
|-----|---------|------------------|
| (30) | §4.1 | Σ(a_ini) = ρ_X(a_ini) / (3 c_D) — initial condition |
| (31) | §4.1 | p_X = w_X · ρ_X,  w_X = −1 + 1/(3 c_D) — pressure |
| (32) | §4.2 | r = ρ_c / ρ_X — density ratio |
| (35) | §4.2 | Ψ(r) = 4 r (r−1) / (1+r)³ — signed kernel |
| (38) | §4.3 | dΣ / dloga = −(1/c_D)[1 − β Ψ(r)] Σ |
| (43) | §4.3 | dρ_X / dloga = −(1/c_D)[1 − β Ψ(r)] ρ_X |
| (55) | §5.1 | ρ_X ∝ a^(−1/c_D) — initial condition on the non-crossing branch |

## Approximation carried forward to Phase 1b.2  →  NOW ENABLED

The CDM back-reaction term in Eq. 44,

    dρ_c / dloga = −3 ρ_c − (β/c_D) Ψ(r) ρ_X,

is wired in by **Phase 1b.2** — see `PHASE1B2_PATCHES.md` (Patch 7
adds a one-line `afterglow_glue_cdm_source()` call inside
`background_derivs()`). Bianchi conservation of the total dark
sector is pinned by Test 11 at 1e-14.

---

## Files added in Phase 1b (no CLASS core edits yet)

    include/afterglow/afterglow_class_glue.h
    source/afterglow/afterglow_class_glue.c
    test/test_afterglow_bg.py          (extended with Tests 8 & 9)

These compile cleanly against the existing `include/afterglow/afterglow.h`
and Phase 1a `afterglow.c`; the new file only depends on `math.h` and
the existing Phase 1a symbols.

---

## Patches to CLASS core files

Five small, surgical patches. Each one is a CLASS-idiomatic
`has_fld`-style insertion and adds no more than ~8 lines.

### Patch 1 — `include/background.h`, top of struct background (near line 112)

Add after the `Omega0_fld` / `use_ppf` / `w0_fld` block:

```c
  /* ─── Afterglow dark energy (Martin & Koh 2026) ─────────────── */
  short has_afterglow;           /**< flag: enable hidden-YM afterglow module */
  double Omega0_X;               /**< afterglow fluid density today           */
  struct afterglow_params ag;    /**< {c_D, beta_aft, Sigma_today, ...}       */
```

Plus add these new integer index fields next to `index_bg_rho_fld` (near line 169):

```c
  int index_bg_rho_X;            /**< afterglow fluid density          */
  int index_bg_p_X;              /**< afterglow fluid pressure         */
  int index_bg_Sigma;            /**< MIS memory variable              */
```

And next to `index_bi_rho_fld` (near line 257):

```c
  int index_bi_rho_X;            /**< {B} afterglow rho_X              */
  int index_bi_Sigma;            /**< {B} afterglow Sigma              */
```

Finally at the top of the file (with the other `#include`s):

```c
#include "afterglow/afterglow.h"
#include "afterglow/afterglow_class_glue.h"
```

### Patch 2 — `source/background.c`, `background_indices()` (~line 963)

Right after the existing `class_define_index(pba->index_bg_rho_fld, pba->has_fld, ...)`
block, add:

```c
  class_define_index(pba->index_bg_rho_X, pba->has_afterglow, index_bg, 1);
  class_define_index(pba->index_bg_p_X,   pba->has_afterglow, index_bg, 1);
  class_define_index(pba->index_bg_Sigma, pba->has_afterglow, index_bg, 1);
```

And next to `class_define_index(pba->index_bi_rho_fld, ...)` (~line 1166):

```c
  class_define_index(pba->index_bi_rho_X, pba->has_afterglow, index_bi, 1);
  class_define_index(pba->index_bi_Sigma, pba->has_afterglow, index_bi, 1);
```

### Patch 3 — `source/background.c`, `background_functions()` (~line 555)

Directly after the `if (pba->has_fld == _TRUE_) { ... }` block that adds
`rho_fld` to `rho_tot`, insert:

```c
  if (pba->has_afterglow == _TRUE_) {
    pvecback[pba->index_bg_rho_X] = pvecback_B[pba->index_bi_rho_X];
    pvecback[pba->index_bg_Sigma] = pvecback_B[pba->index_bi_Sigma];
    double p_X;
    afterglow_glue_pressure(&(pba->ag),
                            pvecback[pba->index_bg_rho_X],
                            &p_X);
    pvecback[pba->index_bg_p_X] = p_X;
    rho_tot += pvecback[pba->index_bg_rho_X];
    p_tot   += p_X;
  }
```

### Patch 4 — `source/background.c`, `background_derivs()` (~line 2660)

Directly after the `if (pba->has_fld == _TRUE_) { dy[pba->index_bi_rho_fld] = ... }`
block, insert:

```c
  if (pba->has_afterglow == _TRUE_) {
    double dX, dS;
    afterglow_glue_derivs(&(pba->ag),
                          pvecback[pba->index_bg_rho_cdm],  /* rho_c */
                          y[pba->index_bi_rho_X],
                          y[pba->index_bi_Sigma],
                          &dX, &dS);
    dy[pba->index_bi_rho_X] = dX;
    dy[pba->index_bi_Sigma] = dS;
  }
```

### Patch 5 — `source/background.c`, `background_solve()` initial conditions (~line 2255)

Directly after the line that sets `pvecback_integration[pba->index_bi_rho_fld] = ...`,
insert:

```c
  if (pba->has_afterglow == _TRUE_) {
    double rX_ini, S_ini;
    afterglow_glue_initial_conditions(&(pba->ag),
                                      pba->Omega0_X,
                                      exp(loga_ini),        /* a_ini */
                                      &rX_ini, &S_ini);
    pvecback_integration[pba->index_bi_rho_X] = rX_ini;
    pvecback_integration[pba->index_bi_Sigma] = S_ini;
  }
```

### Patch 6 — `source/input.c`, parameter parsing

In `input_read_parameters_species()` (search for the `fluid_equation_of_state`
block), add a new block:

```c
  /* ─── Afterglow dark energy (Martin & Koh 2026) ─────────────── */
  class_call(parser_read_double(pfc, "afterglow_on", &param1, &flag1, errmsg),
             errmsg, errmsg);
  if (flag1 == _TRUE_ && param1 > 0.) {
    pba->has_afterglow = _TRUE_;
    pba->ag.afterglow_on = 1;

    /* default parameter values */
    pba->Omega0_X        = 0.69;
    pba->ag.c_D          = 1.0;
    pba->ag.beta_aft     = 0.0;
    pba->ag.Sigma_today  = -1.0;     /* -1 => auto: Omega0_X/(3 c_D) */
    pba->ag.cs2_X        = 1.0;
    pba->ag.z_conf       = 1.0e5;
    pba->ag.Omega0_dr_dark = 0.0;

    class_call(parser_read_double(pfc, "Omega_X",   &param1, &flag1, errmsg), errmsg, errmsg);
    if (flag1 == _TRUE_) pba->Omega0_X     = param1;

    class_call(parser_read_double(pfc, "c_D",       &param1, &flag1, errmsg), errmsg, errmsg);
    if (flag1 == _TRUE_) pba->ag.c_D       = param1;

    class_call(parser_read_double(pfc, "beta_aft",  &param1, &flag1, errmsg), errmsg, errmsg);
    if (flag1 == _TRUE_) pba->ag.beta_aft  = param1;

    class_call(parser_read_double(pfc, "Sigma_today", &param1, &flag1, errmsg), errmsg, errmsg);
    if (flag1 == _TRUE_) pba->ag.Sigma_today = param1;

    /* disable CLASS's default lambda so total Omega sums to 1 */
    pba->Omega0_lambda = 0.;
  }
  else {
    pba->has_afterglow = _FALSE_;
  }
```

---

## Verification

After applying the six patches and running `make class`, the following
must all hold:

1. **LCDM regression.** `./class explanatory.ini` (with no `afterglow_on`
   line) must produce bitwise-identical background table to the unmodified
   CLASS. This is guaranteed by the `has_afterglow == _FALSE_` gate on every
   new block.

2. **β = 0 regression.** `./class afterglow.ini` with `beta_aft = 0` and
   `c_D` very large should recover LCDM background densities to 1e-10
   relative.

3. **Unit tests.** `python test/test_afterglow_bg.py` must report
   **24 passed, 0 failed** (Tests 8 and 9 specifically pin the glue layer
   point-wise against the Phase 1a rhs).

4. **Phase 1a / 1b agreement.** The rho_X table returned by CLASS's own
   evolver (after Patch 4) must agree with `afterglow_background_evolve()`
   to the CLASS background tolerance (~1e-8) on the non-crossing benchmark.

## Rollback

Every Phase 1b change is gated by `if (pba->has_afterglow == _TRUE_)`. To
disable the module entirely, simply omit `afterglow_on` from the `.ini`
file or set it to zero. No recompile needed.

## Next phase (1b.2)

Phase 1b.2 adds the CDM back-reaction term from Eq. 44 directly to
`dy[pba->index_bi_rho_cdm]` in Patch 4, and re-runs both the unit tests
and the CLASS LCDM regression with `beta_aft = 0` to confirm nothing
moves on the β = 0 branch.

## Next phase (2)

Phase 2 adds perturbations: δ_X, θ_X, δΣ with the causal sound speed
c²_X = 1 and the MIS memory variable δΣ as an **independent** degree of
freedom (do not fold into δ_X). The perturbation equations are Eqs. (?)
in §6 of the paper; see `PLAN_AFTERGLOW.md` for the rollout.
