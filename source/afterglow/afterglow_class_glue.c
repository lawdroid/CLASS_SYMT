/** @file afterglow_class_glue.c
 *
 *  Phase 1b glue between the self-contained afterglow module (Phase 1a)
 *  and CLASS's background.c / input.c. See afterglow_class_glue.h for
 *  the full design note and equation references.
 *
 *  This file contains ONLY plain C helpers that implement paper
 *  Eqs. 38, 43, 31, 30, and 55 (all VERIFIED). No CLASS state vector
 *  layout is referenced here — background.c is the only place where
 *  pvecback_B / pvecback indices live.
 *
 *  Paper:  Martin & Koh (April 2026),
 *          "Dark Energy as the Thermodynamic Afterglow of a Hidden
 *           Gauge-Theory Transition"  (docs/Martin_Koh_Afterglow_April2026.pdf)
 */

#include <math.h>
#include "afterglow/afterglow.h"
#include "afterglow/afterglow_class_glue.h"

/* ------------------------------------------------------------------ */
/*  Derivative helper — Eqs. 38, 43                                    */
/*                                                                     */
/*   d rho_X / d(log a) = -(1/c_D) [1 - beta Psi(r)] rho_X   (Eq. 43)  */
/*   d Sigma / d(log a) = -(1/c_D) [1 - beta Psi(r)] Sigma   (Eq. 38)  */
/*                                                                     */
/*   r = rho_c / rho_X   (Eq. 32)                                      */
/*   Psi(r) = 4 r (r - 1) / (1 + r)^3   (Eq. 35)                       */
/* ------------------------------------------------------------------ */
int afterglow_glue_derivs(struct afterglow_params * ap,
                          double rho_c,
                          double rho_X, double Sigma,
                          double * d_rho_X, double * d_Sigma) {
  double r;
  double Psi;
  double factor;

  /* r = rho_c / rho_X, guarded */
  r = (rho_X > 1e-300) ? (rho_c / rho_X) : 0.0;

  /* Psi(r) — uses the Phase 1a kernel */
  afterglow_kernel_Psi(r, &Psi);

  /* Common (1/c_D) [1 - beta * Psi(r)] prefactor */
  factor = (1.0 / ap->c_D) * (1.0 - ap->beta_aft * Psi);

  *d_rho_X = -factor * rho_X;   /* Eq. 43 */
  *d_Sigma = -factor * Sigma;   /* Eq. 38 */

  return _SUCCESS_;
}

/* ------------------------------------------------------------------ */
/*  Pressure helper — Eq. 31                                           */
/*                                                                     */
/*   w_X = -1 + 1/(3 c_D)                                              */
/*   p_X = w_X * rho_X                                                 */
/* ------------------------------------------------------------------ */
int afterglow_glue_pressure(struct afterglow_params * ap,
                            double rho_X, double * p_X) {
  double w_X = -1.0 + 1.0 / (3.0 * ap->c_D);
  *p_X = w_X * rho_X;
  return _SUCCESS_;
}

/* ------------------------------------------------------------------ */
/*  Initial conditions at the earliest time — Eq. 55 + Eq. 30          */
/*                                                                     */
/*   rho_X(a_ini) ≈ Omega0_X * a_ini^(-1/c_D)     (Eq. 55, beta=0)     */
/*   Sigma(a_ini) = rho_X(a_ini) / (3 c_D)        (Eq. 30)             */
/*                                                                     */
/*  This is exact in the beta=0 limit and an excellent starting guess  */
/*  on the non-crossing branch. Because |beta * Psi(r)| ≤ 2 sqrt(3)/9  */
/*  |beta| ~ 0.08 |beta|, the leading backward error is of order       */
/*  0.1 |beta| e-folds from the crossing region.                       */
/* ------------------------------------------------------------------ */
int afterglow_glue_initial_conditions(struct afterglow_params * ap,
                                      double Omega0_X, double a_ini,
                                      double * rho_X_ini,
                                      double * Sigma_ini) {
  /* log(a_ini) < 0; a_ini^(-1/c_D) = exp(-log(a_ini)/c_D) */
  double loga_ini = log(a_ini);
  *rho_X_ini = Omega0_X * exp(-loga_ini / ap->c_D);
  *Sigma_ini = (*rho_X_ini) / (3.0 * ap->c_D);
  return _SUCCESS_;
}

/* ------------------------------------------------------------------ */
/*  Phase 1b.2 — CDM back-reaction source (Eq. 44)                    */
/*                                                                     */
/*   d rho_c / d(log a) = -3 rho_c  +  cdm_source                     */
/*                                                                     */
/*   cdm_source = -(beta / c_D) Psi(r) rho_X                          */
/*                                                                     */
/*  Bianchi check (total dark sector):                                */
/*     d/dN (rho_X + rho_c) + 3 rho_c + 3 (1+w_X) rho_X               */
/*     = +(beta/c_D) Psi rho_X  +  [-(beta/c_D) Psi rho_X]            */
/*     = 0                                                            */
/*  so the pair (Eq. 43, Eq. 44) conserves total energy modulo the    */
/*  standard Hubble-drag terms.                                       */
/* ------------------------------------------------------------------ */
int afterglow_glue_cdm_source(struct afterglow_params * ap,
                              double rho_c, double rho_X,
                              double * cdm_source) {
  double r;
  double Psi;

  /* r = rho_c / rho_X, guarded */
  r = (rho_X > 1e-300) ? (rho_c / rho_X) : 0.0;

  /* Psi(r) — uses the Phase 1a kernel */
  afterglow_kernel_Psi(r, &Psi);

  /* Eq. 44 source term (beyond -3 rho_c) */
  *cdm_source = -(ap->beta_aft / ap->c_D) * Psi * rho_X;

  return _SUCCESS_;
}
