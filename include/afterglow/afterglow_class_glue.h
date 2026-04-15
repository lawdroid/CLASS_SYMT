/** @file afterglow_class_glue.h
 *
 *  Phase 1b glue layer between the (self-contained) afterglow module
 *  of Phase 1a and CLASS's background.c / input.c.
 *
 *  The goal of this file is to keep all new physics out of background.c.
 *  The only pieces that background.c sees are four function-call
 *  one-liners — one in background_indices, one in background_functions,
 *  one in background_derivs, one in background_solve (initial
 *  conditions) — plus four new index fields in `struct background`.
 *
 *  Paper equations used (VERIFIED, Martin & Koh, April 2026):
 *
 *      (31)  w_X = -1 + 1/(3 c_D)                                 [§4.1]
 *      (32)  r   = rho_c / rho_X                                   [§4.2]
 *      (35)  Psi(r) = 4 r (r-1) / (1+r)^3                          [§4.2]
 *      (43)  d rho_X / d(log a) = -(1/c_D)[1 - beta Psi(r)] rho_X [§4.3]
 *      (38)  d Sigma / d(log a) = -(1/c_D)[1 - beta Psi(r)] Sigma [§4.3]
 *      (30)  Sigma_today = rho_X_today / (3 c_D)                  [§4.1]
 *
 *  Phase 1b.2 ADDITION:
 *      (44)  d rho_c / d(log a) = -3 rho_c - (beta/c_D) Psi(r) rho_X
 *      The CDM back-reaction source term is now exposed via
 *      afterglow_glue_cdm_source() so background.c can add it to
 *      dy[pba->index_bi_rho_cdm]. The source is written as the
 *      EXTRA term beyond the standard -3 rho_c, so the call site
 *      stays:
 *          dy[pba->index_bi_rho_cdm] += source;
 *      with source = -(beta/c_D) Psi(r) rho_X  (Eq. 44).
 *
 *      Bianchi / total-dark-sector conservation pins this sign:
 *      the matching source to rho_X is +(beta/c_D) Psi(r) rho_X
 *      (already inside Eq. 43 via the [1 - beta Psi] factor), so
 *      d(rho_X + rho_c)/d(log a) is source-free.
 */

#ifndef __AFTERGLOW_CLASS_GLUE__
#define __AFTERGLOW_CLASS_GLUE__

#include "afterglow/afterglow.h"

/* Forward decl so we don't drag all of background.h into this header. */
struct background;

/**
 *  Compute the derivative of (rho_X, Sigma) with respect to log(a),
 *  given the current log(a), the current CDM density, and the current
 *  values of rho_X and Sigma on the background integration grid.
 *
 *  Implements Eqs. 38 and 43 point-wise. Safe to call from inside
 *  background_derivs() without knowing anything about the CLASS
 *  background state vector layout.
 *
 *  @param ap       [in]  afterglow parameter block
 *  @param rho_c    [in]  current CDM energy density (same units as rho_X)
 *  @param rho_X    [in]  current rho_X from y[index_bi_rho_X]
 *  @param Sigma    [in]  current Sigma from y[index_bi_Sigma]
 *  @param d_rho_X  [out] d rho_X / d(log a)   (Eq. 43)
 *  @param d_Sigma  [out] d Sigma / d(log a)   (Eq. 38)
 */
int afterglow_glue_derivs(struct afterglow_params * ap,
                          double rho_c,
                          double rho_X, double Sigma,
                          double * d_rho_X, double * d_Sigma);

/**
 *  Compute the pressure p_X from rho_X via the intrinsic EoS (Eq. 31).
 *  Pressure feeds directly into CLASS's p_tot bookkeeping.
 *
 *  @param ap     [in]  afterglow parameter block
 *  @param rho_X  [in]  current afterglow fluid density
 *  @param p_X    [out] pressure such that p_X = w_X * rho_X, Eq. 31
 */
int afterglow_glue_pressure(struct afterglow_params * ap,
                            double rho_X, double * p_X);

/**
 *  Compute initial conditions at the *earliest* time in the CLASS
 *  background integration (large negative log(a)). CLASS integrates
 *  forward in log(a) so initial conditions must be supplied at the
 *  LEFT edge of the grid.
 *
 *  We are given the target today values (Omega0_X at a=1) and must
 *  evolve Eq. 43 backward analytically or assume a simple scaling to
 *  supply (rho_X_ini, Sigma_ini) at a_ini. In Phase 1b we use the
 *  beta=0 analytic backward extrapolation:
 *
 *      rho_X_ini = Omega0_X * a_ini^(-1/c_D)     (Eq. 55)
 *      Sigma_ini = rho_X_ini / (3 c_D)           (Eq. 30)
 *
 *  This is exact in the beta=0 limit and is the natural Phase 1b
 *  starting guess. Phase 1b.2 will refine this with a shooting step
 *  if needed.
 *
 *  @param ap          [in]  afterglow parameter block
 *  @param Omega0_X    [in]  rho_X today in units of rho_crit_0
 *  @param a_ini       [in]  scale factor at start of CLASS integration
 *  @param rho_X_ini   [out] initial rho_X at a = a_ini
 *  @param Sigma_ini   [out] initial Sigma at a = a_ini
 */
int afterglow_glue_initial_conditions(struct afterglow_params * ap,
                                      double Omega0_X, double a_ini,
                                      double * rho_X_ini,
                                      double * Sigma_ini);

/**
 *  Phase 1b.2 — CDM back-reaction source term from Eq. 44.
 *
 *      d rho_c / d(log a) = -3 rho_c - (beta/c_D) Psi(r) rho_X
 *
 *  This helper returns ONLY the EXTRA term beyond the standard
 *  -3 rho_c pressureless-matter derivative. The CLASS call site
 *  stays:
 *
 *      dy[pba->index_bi_rho_cdm] += cdm_source;   (Eq. 44)
 *
 *  Gated by has_afterglow in background.c so the LCDM regression
 *  path is bit-identical.
 *
 *  @param ap          [in]  afterglow parameter block
 *  @param rho_c       [in]  current CDM density (for r = rho_c/rho_X)
 *  @param rho_X       [in]  current afterglow fluid density
 *  @param cdm_source  [out] -(beta/c_D) Psi(r) rho_X     (Eq. 44)
 */
int afterglow_glue_cdm_source(struct afterglow_params * ap,
                              double rho_c, double rho_X,
                              double * cdm_source);

#endif /* __AFTERGLOW_CLASS_GLUE__ */
