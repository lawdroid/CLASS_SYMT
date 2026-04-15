# Phase 2 — Afterglow linear perturbations (δ_X, θ_X, σ̂) with MIS closure

Phase 1 (background) and 1b + 1b.2 (CLASS background wiring + CDM
back-reaction) are complete. Phase 2 adds three independent
perturbation DOFs for the afterglow fluid onto the synchronous-gauge
Boltzmann grid — kept OUT of CLASS's stock `fld` module, which
collapses the memory variable Σ into a barotropic closure and is
unsuitable for Maxwell-Israel-Stewart bulk stress.

All new physics lives in:

    include/afterglow/afterglow_pert.h
    source/afterglow/afterglow_pert.c
    test/test_afterglow_pt.py

`perturbations.c` only ever calls three glue helpers. No new
perturbation physics lives in CLASS core files.

---

## Variables

| Symbol     | Meaning                                               |
|------------|-------------------------------------------------------|
| `δ_X`      | δρ_X / ρ̄_X (density contrast)                        |
| `θ_X`      | k_i δu^i / a (velocity divergence, CLASS convention)  |
| `σ̂`       | δΣ / Σ̄ (fractional MIS memory perturbation)          |

Pressure closure:

    δp_X = w_X · ρ̄_X · σ̂          (MIS — ties δp to δΣ, NOT δρ)

so the effective sound speed δp_X / δρ_X = w_X · σ̂/δ_X is **not**
locked to w_X. It collapses to w_X only on the adiabatic sub-manifold
σ̂ = δ_X.

---

## Evolution equations (synchronous gauge, conformal time)

With `r̄ = ρ̄_c/ρ̄_X`, `Ψ(r) = 4r(r-1)/(1+r)³`,
`Ψ'(r) = 4(-r²+4r-1)/(1+r)⁴`, `w_X = -1 + 1/(3c_D)`,
`δΨ = Ψ'(r̄) · r̄ · (δ_c - δ_X)`:

    (P1)  δ_X' = -(1+w_X)(θ_X + h'/2)
                 - 3𝓗 w_X (σ̂ - δ_X)                ← MIS bulk
                 + 𝓗 (β/c_D) δΨ                    ← δQ from δΨ

    (P2)  θ_X' = -𝓗 (1 - 3w_X) θ_X
                 + (k² w_X / (1+w_X)) · σ̂          ← pressure gradient

    (P3)  σ̂'  = (𝓗/c_D) β r̄ Ψ'(r̄) (δ_c - δ_X)
                 - ((θ_X + h'/2) / (3c_D)) (1 - β Ψ̄)

In the LCDM limit (β=0, c_D→∞): 1+w_X → 0, all drivers decouple, and
the afterglow fluid plays no role in perturbations.

---

## Patches to CLASS core files

### Patch 8 — `include/perturbations.h` (near fld indices, ~line 250)

```c
  /* ─── Afterglow perturbations (Martin & Koh 2026) ──────── */
  int index_pt_delta_X;      /**< afterglow density contrast    */
  int index_pt_theta_X;      /**< afterglow velocity divergence */
  int index_pt_sigma_hat;    /**< afterglow MIS memory perturbation */
```

### Patch 9 — `perturb_indices_of_current_vectors()` (next to fld slots)

```c
  if (pba->has_afterglow == _TRUE_) {
    class_define_index(ppw->pv->index_pt_delta_X,    _TRUE_, index_pt, 1);
    class_define_index(ppw->pv->index_pt_theta_X,    _TRUE_, index_pt, 1);
    class_define_index(ppw->pv->index_pt_sigma_hat,  _TRUE_, index_pt, 1);
  }
```

### Patch 10 — `perturb_initial_conditions()` (adiabatic mode block)

```c
  if (pba->has_afterglow == _TRUE_) {
    double dX, tX, sH;
    class_call(afterglow_pert_initial_conditions(&(pba->ag),
                                                 ppw->pvecmetric[ppw->index_mt_delta_g]
                                                   /* δ_γ at IC */,
                                                 ppw->pvecmetric[ppw->index_mt_theta_g]
                                                   /* θ_γ at IC */,
                                                 &dX, &tX, &sH),
               errmsg, errmsg);
    y[ppw->pv->index_pt_delta_X]   = dX;
    y[ppw->pv->index_pt_theta_X]   = tX;
    y[ppw->pv->index_pt_sigma_hat] = sH;
  }
```

### Patch 11 — `perturb_derivs()` (new block after fld)

```c
  if (pba->has_afterglow == _TRUE_) {
    double d_delta_X, d_theta_X, d_sigma_hat;
    class_call(afterglow_pert_rhs(&(pba->ag),
                                  pvecback[pba->index_bg_H_conf]  /* aH */,
                                  k,
                                  pvecback[pba->index_bg_rho_cdm],
                                  y[pba->index_bi_rho_X],
                                  y[ppw->pv->index_pt_delta_cdm],
                                  y[ppw->pv->index_pt_delta_X],
                                  y[ppw->pv->index_pt_theta_X],
                                  y[ppw->pv->index_pt_sigma_hat],
                                  ppw->pvecmetric[ppw->index_mt_h_prime],
                                  &d_delta_X, &d_theta_X, &d_sigma_hat),
               errmsg, errmsg);
    dy[ppw->pv->index_pt_delta_X]   = d_delta_X;
    dy[ppw->pv->index_pt_theta_X]   = d_theta_X;
    dy[ppw->pv->index_pt_sigma_hat] = d_sigma_hat;
  }
```

### Patch 12 — `perturb_einstein()` (source of δp_tot and δρ_tot)

```c
  if (pba->has_afterglow == _TRUE_) {
    double delta_p_X;
    class_call(afterglow_pert_pressure(&(pba->ag),
                                       pvecback[pba->index_bg_rho_X],
                                       y[ppw->pv->index_pt_sigma_hat],
                                       &delta_p_X),
               errmsg, errmsg);
    ppw->delta_rho_tot += y[pba->index_bi_rho_X] * y[ppw->pv->index_pt_delta_X];
    ppw->delta_p_tot   += delta_p_X;
    ppw->rho_plus_p_theta_tot +=
        (1.0 + (-1.0 + 1.0/(3.0 * pba->ag.c_D)))
          * pvecback[pba->index_bg_rho_X]
          * y[ppw->pv->index_pt_theta_X];
  }
```

All patches are gated by `has_afterglow == _TRUE_` so the LCDM
regression stays bit-identical.

---

## Regression coverage (`python test/test_afterglow_pt.py`)

    Test 12  LCDM limit (c_D→∞, β=0) — all RHS finite & expected
    Test 13  Adiabatic sub-manifold fixed point (all RHS = 0)
    Test 14  Psi'(r) matches numerical derivative of Psi(r) to 1e-8
    Test 15  MIS pressure closure and off-sub-manifold deviation
    Test 16  Memory relaxation coefficient = 1/(3 c_D) per e-fold
    Test 17  Psi' coupling coefficient match between (P1) and (P3)
    Test 18  Adiabatic ICs: δ_X = δ_γ/(4 c_D), σ̂ = δ_X, θ_X = θ_γ

    RESULT:  27 passed, 0 failed

Phase 1 + 1b + 1b.2 regression (`test/test_afterglow_bg.py`):
28 passed, 0 failed.

---

## What still defers to Phase 3

* **MCMC** over `{z_c, ΔN_eff, c_D}` + standard ΛCDM on Planck 2018 +
  DESI DR2 + Pantheon+.
* **Shooting** for Ω_X today (currently beta=0 analytic extrapolation
  from a_ini back to today).
* **Momentum exchange perturbation** (δQ_momentum ≡ 0 in Phase 2 — the
  4-velocity frame is locked to u_c, so no linear-order momentum
  transfer to X).
* **Anisotropic stress** (σ^{aniso}_X ≡ 0 in Phase 2 — only the BULK
  MIS stress is active; shear memory would need a separate relaxation
  equation, not present in Martin & Koh's derivation).
