/** @file afterglow_pert.h
 *
 *  Phase 2 — Linear perturbation glue for the afterglow dark-energy
 *  fluid of Martin & Koh (April 2026). Implements three new DOFs on
 *  the synchronous-gauge Boltzmann grid:
 *
 *      delta_X  = delta rho_X / rho_X_bg       (density contrast)
 *      theta_X                                 (velocity divergence)
 *      sigma_hat = delta Sigma / Sigma_bg      (MIS memory perturbation)
 *
 *  Kept STRICTLY OUT of CLASS's stock `fld` module because the fld
 *  module collapses the memory variable and thus assigns a barotropic
 *  closure `delta p = c_s^2 delta rho`. The afterglow closure is
 *  instead
 *
 *      delta p_X = w_X * rho_X_bg * sigma_hat
 *
 *  which is tied to delta_Sigma independently of delta_X, exactly as
 *  the paper's Maxwell-Israel-Stewart bulk-stress memory requires.
 *
 *  Paper equations (§4.3, §4.4) VERIFIED against
 *  docs/Martin_Koh_Afterglow_April2026.pdf :
 *
 *      (31)  w_X       = -1 + 1/(3 c_D)                          [§4.1]
 *      (32)  r         = rho_c / rho_X                            [§4.2]
 *      (35)  Psi(r)    = 4 r (r-1) / (1+r)^3                      [§4.2]
 *      (37)  u^mu d_mu Sigma = -(Theta/(3 c_D))[1 - beta Psi] Sigma   [§4.3]
 *      (42)  Q         = -(beta/c_D) H Psi(r) rho_X               [§4.3]
 *
 *  Phase 2 perturbation equations (synchronous gauge, conformal time):
 *
 *      (P1)  delta_X'  = -(1+w_X)(theta_X + h'/2)
 *                        - 3 H_conf w_X (sigma_hat - delta_X)
 *                        + H_conf (beta/c_D) delta_Psi
 *
 *      (P2)  theta_X'  = -H_conf (1 - 3 w_X) theta_X
 *                        + (k^2 w_X / (1 + w_X)) sigma_hat
 *
 *      (P3)  sigma_hat' = (H_conf / c_D) beta r_bar Psi_prime(r_bar)
 *                                                  (delta_c - delta_X)
 *                         - ((theta_X + h'/2) / (3 c_D))(1 - beta Psi_bar)
 *
 *  Auxiliary / closure:
 *
 *      delta_Psi    = Psi_prime(r_bar) * r_bar * (delta_c - delta_X)
 *      Psi_prime(r) = 4 (-r^2 + 4 r - 1) / (1 + r)^4
 *      delta p_X    = w_X * rho_X_bg * sigma_hat        (MIS closure)
 *
 *  In the LCDM limit (c_D -> inf, beta=0) all three RHS vanish, so
 *  the afterglow perturbations decouple (Test 12).
 */

#ifndef __AFTERGLOW_PERT__
#define __AFTERGLOW_PERT__

#include "afterglow/afterglow.h"

/* Forward decls so we avoid pulling CLASS headers into this file. */
struct afterglow_params;

/* Indices for the three new perturbation DOFs. Must match the order
   in which background.c / perturbations.c reserves slots via
   class_define_index for this fluid. */
enum afterglow_pert_idx {
  AFTERGLOW_PT_DELTA_X = 0,   /* delta rho_X / rho_X_bg */
  AFTERGLOW_PT_THETA_X = 1,   /* velocity divergence    */
  AFTERGLOW_PT_SIGMA_H = 2,   /* delta Sigma / Sigma_bg */
  AFTERGLOW_PT_N       = 3
};

/**
 *  Psi'(r) — analytic derivative of the signed loading kernel.
 *  Appears in delta_Psi via  delta_Psi = Psi'(r_bar) * delta_r,
 *  with delta_r = r_bar * (delta_c - delta_X).
 *
 *      Psi'(r) = 4 (-r^2 + 4 r - 1) / (1 + r)^4
 *
 *  Zeros of Psi'(r) are the extrema of Psi(r), at r = 2 +- sqrt(3),
 *  where |Psi| = 2 sqrt(3)/9.
 */
int afterglow_kernel_Psi_prime(double r, double * Psi_prime);

/**
 *  Perturbation RHS packed into a single call. All inputs are
 *  BACKGROUND quantities at the current conformal time plus the
 *  PERTURBATION state vector (delta_X, theta_X, sigma_hat). Outputs
 *  the three time derivatives in synchronous gauge, conformal time.
 *
 *  All new physics lives here; perturbations.c only calls this glue.
 *
 *  @param ap           [in]  afterglow parameter block
 *  @param H_conf       [in]  conformal Hubble = a*H
 *  @param k            [in]  comoving wavenumber (1/Mpc)
 *  @param rho_c_bg     [in]  background rho_c at this conformal time
 *  @param rho_X_bg     [in]  background rho_X at this conformal time
 *  @param delta_c      [in]  CDM density contrast (for coupling)
 *  @param delta_X      [in]  afterglow density contrast
 *  @param theta_X      [in]  afterglow velocity divergence
 *  @param sigma_hat    [in]  afterglow MIS memory perturbation
 *  @param h_prime_syn  [in]  h' in synchronous-gauge convention
 *  @param d_delta_X    [out] (P1)
 *  @param d_theta_X    [out] (P2)
 *  @param d_sigma_hat  [out] (P3)
 */
int afterglow_pert_rhs(struct afterglow_params * ap,
                       double H_conf, double k,
                       double rho_c_bg, double rho_X_bg,
                       double delta_c,
                       double delta_X, double theta_X, double sigma_hat,
                       double h_prime_syn,
                       double * d_delta_X,
                       double * d_theta_X,
                       double * d_sigma_hat);

/**
 *  MIS pressure-perturbation closure (for perturbations.c to feed
 *  into the metric Einstein equations):
 *
 *      delta p_X = w_X * rho_X_bg * sigma_hat                   (P-closure)
 *
 *  Note the explicit tie to sigma_hat (NOT delta_X). This is the
 *  signature of the MIS bulk-stress memory and the reason the stock
 *  CLASS fld module is unsuitable.
 */
int afterglow_pert_pressure(struct afterglow_params * ap,
                            double rho_X_bg, double sigma_hat,
                            double * delta_p_X);

/**
 *  Synchronous-gauge adiabatic initial conditions for the three
 *  perturbation DOFs, far inside the radiation era (where the
 *  adiabatic mode dominates). The standard CLASS prescription
 *  gives all species the adiabatic ratio:
 *
 *      delta_i / (1 + w_i) = delta_r / (1 + w_r) = ... = -h/2 ... etc.
 *
 *  For the afterglow fluid the memory variable sigma_hat is
 *  initialized to the equilibrium value sigma_hat = delta_X so the
 *  MIS closure starts in the barotropic sub-manifold (pressure
 *  perturbation matches the adiabatic sound speed). Any departure
 *  from sigma_hat = delta_X at later times is then the DIAGNOSTIC
 *  signature of the MIS memory.
 *
 *  theta_X is initialized to the adiabatic velocity at IC time.
 *
 *  @param ap             [in]  afterglow parameter block
 *  @param delta_gamma_IC [in]  photon density contrast at IC time
 *  @param theta_gamma_IC [in]  photon velocity divergence at IC time
 *  @param delta_X_ini    [out] IC for delta_X    (adiabatic)
 *  @param theta_X_ini    [out] IC for theta_X    (adiabatic)
 *  @param sigma_hat_ini  [out] IC for sigma_hat  (= delta_X_ini)
 */
int afterglow_pert_initial_conditions(struct afterglow_params * ap,
                                      double delta_gamma_IC,
                                      double theta_gamma_IC,
                                      double * delta_X_ini,
                                      double * theta_X_ini,
                                      double * sigma_hat_ini);

#endif /* __AFTERGLOW_PERT__ */
