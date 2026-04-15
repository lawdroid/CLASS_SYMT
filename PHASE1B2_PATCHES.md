# Phase 1b.2 — CDM back-reaction (Eq. 44) wired into `background_derivs()`

Phase 1b wired the afterglow fluid `(rho_X, Sigma)` into CLASS's
`background.c` via four glue calls, with the CDM back-reaction term of
Eq. 44 deferred as a documented approximation. Phase 1b.2 now enables
that term.

Paper (all equation numbers verified against
`docs/Martin_Koh_Afterglow_April2026.pdf`):

| Eq. | Section | Role in Phase 1b.2 |
|-----|---------|--------------------|
| (42) | §4.3 | `Q = -(β/c_D) H Ψ(r) ρ_X` — exchange current |
| (43) | §4.3 | `dρ_X / dN = -(1/c_D)[1 - β Ψ] ρ_X` (already in Phase 1b) |
| (44) | §4.3 | `dρ_c / dN = -3 ρ_c - (β/c_D) Ψ(r) ρ_X`   **← added here** |

Bianchi identity (total dark sector):

    d/dN (ρ_X + ρ_c) + 3 ρ_c + 3 (1 + w_X) ρ_X = 0

with `(1 + w_X) = 1/(3 c_D)`. Tested at 1e-14 by Test 11 in
`test/test_afterglow_bg.py`.

---

## Files changed in Phase 1b.2

Additive only — no existing line is edited:

    include/afterglow/afterglow_class_glue.h
        + prototype: afterglow_glue_cdm_source()
        - doc block: "Phase 1b.2 APPROXIMATION" note replaced by
                     "Phase 1b.2 ADDITION" description of Eq. 44

    source/afterglow/afterglow_class_glue.c
        + implementation: afterglow_glue_cdm_source()
          returns the EXTRA term beyond -3 ρ_c:
              source = -(β/c_D) Ψ(r) ρ_X

    test/test_afterglow_bg.py
        + glue_cdm_source() Python mirror
        + test10  : point-match to Eq. 44 (rel err 0)
        + test11  : total-dark-sector Bianchi conservation (rel err 8e-17)

Phase 1a (`source/afterglow/afterglow.c`) is **unchanged** — the self-
contained integrator keeps the Phase-1 CDM background
`ρ_c(N) = Ω₀_c e^{-3N}` so Tests 1–7 remain tautological regression
guards. Phase 1b.2 lives entirely in the glue layer, which is where
CLASS-side integration will be driven.

---

## Patch 7 — `source/background.c`, inside `background_derivs()`

Add the single extra source term next to the afterglow glue call that
Patch 4 already introduced (Phase 1b). The complete Phase 1b + 1b.2
block reads:

```c
  if (pba->has_afterglow == _TRUE_) {
    double d_rho_X, d_Sigma;
    class_call(afterglow_glue_derivs(&(pba->ag),
                                     pvecback[pba->index_bg_rho_cdm],
                                     y[pba->index_bi_rho_X],
                                     y[pba->index_bi_Sigma],
                                     &d_rho_X, &d_Sigma),
               errmsg, errmsg);
    dy[pba->index_bi_rho_X] = d_rho_X;   /* Eq. 43 */
    dy[pba->index_bi_Sigma] = d_Sigma;   /* Eq. 38 */

    /* ─── Phase 1b.2: CDM back-reaction (Eq. 44) ─────────────── */
    double cdm_source;
    class_call(afterglow_glue_cdm_source(&(pba->ag),
                                         pvecback[pba->index_bg_rho_cdm],
                                         y[pba->index_bi_rho_X],
                                         &cdm_source),
               errmsg, errmsg);
    dy[pba->index_bi_rho_cdm] += cdm_source;   /* Eq. 44 */
  }
```

Nothing else in `background_derivs()` changes. The `+=` is critical:
the standard pressureless-matter derivative `-3 ρ_c` is already
written by CLASS one block above, so we only add the
`-(β/c_D) Ψ(r) ρ_X` correction.

Gated entirely by `has_afterglow == _TRUE_`, so LCDM regression remains
bit-identical.

---

## Expected size of the correction

On the non-crossing benchmark (`β ≲ 1`, matter-dominated past):

    |correction| / (3 ρ_c)  ≈  (β/c_D) · |Ψ|_max · (ρ_X / ρ_c)
                            ≲  0.38 · |β|/c_D · (ρ_X / ρ_c)

At matter–DE equality `r = 1`, `Ψ = 0`, so the correction vanishes.
Peak contribution is at `r = 2 ± √3`, where `|Ψ| = 2√3/9 ≈ 0.385`. In
a typical run with `β = 0.1`, `c_D = 1.3` this is ~3% of the drag
term — the order of magnitude the paper note already anticipated.

---

## Regression coverage

`python test/test_afterglow_bg.py` — 28 passed, 0 failed.

New tests:

    Test 10  cdm_source point-match to Eq. 44          max rel err 0
    Test 11  Bianchi cancellation β/c_D · Ψ · ρ_X      max rel residual 8e-17

Tests 1–9 continue to pass unchanged (24 → 28).

---

## What still defers to later phases

* Phase 2 (perturbations) — `δ_X, θ_X, δΣ` with Maxwell–Cattaneo bulk
  stress, independent of the stock `fld` module.
* Phase 3 — MCMC over `{z_c, ΔN_eff, c_D}` + ΛCDM on Planck + DESI DR2
  + Pantheon+.
