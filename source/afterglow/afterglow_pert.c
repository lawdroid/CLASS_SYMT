/** @file afterglow_pert.c
 *
 *  Phase 2 — Linear perturbation glue for the afterglow fluid
 *  (Martin & Koh, April 2026). See afterglow_pert.h for the full
 *  equation list and design notes.
 *
 *  Units and gauge convention:
 *    * synchronous gauge (CLASS default)
 *    * conformal time tau, H_conf = a H = a' / a
 *    * velocity divergence theta = partial_i v^i in CLASS normalization
 *    * h' is the synchronous-gauge metric trace rate
 *    * k is the comoving wavenumber in 1/Mpc
 *
 *  This file depends ONLY on math.h and on the Phase 1 symbol
 *  afterglow_kernel_Psi() — it does not include any CLASS internals.
 */

#include <math.h>
#include "afterglow/afterglow.h"
#include "afterglow/afterglow_pert.h"

/* ------------------------------------------------------------------ */
/*  Psi'(r) — analytic derivative of Psi(r) = 4 r (r-1) / (1+r)^3     */
/*                                                                     */
/*    Psi'(r) = 4 (-r^2 + 4 r - 1) / (1 + r)^4                        */
/*                                                                     */
/*  Zeros of Psi'(r): r = 2 +- sqrt(3)  (extrema of Psi, |Psi|_max)   */
/* ------------------------------------------------------------------ */
int afterglow_kernel_Psi_prime(double r, double * Psi_prime) {
  double num = 4.0 * (-r * r + 4.0 * r - 1.0);
  double one_plus_r = 1.0 + r;
  double denom = one_plus_r * one_plus_r * one_plus_r * one_plus_r;
  *Psi_prime = num / denom;
  return _SUCCESS_;
}

/* ------------------------------------------------------------------ */
/*  Perturbation RHS — synchronous gauge, conformal time              */
/*                                                                     */
/*  (P1)  delta_X' = -(1+w_X)(theta_X + h'/2)                          */
/*                   - 3 H_conf w_X (sigma_hat - delta_X)              */
/*                   + H_conf (beta/c_D) delta_Psi                     */
/*                                                                     */
/*  (P2)  theta_X' = -H_conf (1 - 3 w_X) theta_X                       */
/*                   + (k^2 w_X / (1 + w_X)) sigma_hat                 */
/*                                                                     */
/*  (P3)  sigma_hat' = (H_conf / c_D) beta r_bar Psi'(r_bar)           */
/*                                        (delta_c - delta_X)          */
/*                    - ((theta_X + h'/2)/(3 c_D)) (1 - beta Psi_bar) */
/*                                                                     */
/*  with delta_Psi = Psi'(r_bar) * r_bar * (delta_c - delta_X),        */
/*       r_bar     = rho_c_bg / rho_X_bg,                              */
/*       w_X       = -1 + 1/(3 c_D).                                   */
/* ------------------------------------------------------------------ */
int afterglow_pert_rhs(struct afterglow_params * ap,
                       double H_conf, double k,
                       double rho_c_bg, double rho_X_bg,
                       double delta_c,
                       double delta_X, double theta_X, double sigma_hat,
                       double h_prime_syn,
                       double * d_delta_X,
                       double * d_theta_X,
                       double * d_sigma_hat) {
  double c_D = ap->c_D;
  double beta = ap->beta_aft;
  double w_X = -1.0 + 1.0 / (3.0 * c_D);
  double one_plus_w = 1.0 + w_X;        /* = 1/(3 c_D),  > 0 strictly */
  double r_bar = (rho_X_bg > 1e-300) ? (rho_c_bg / rho_X_bg) : 0.0;

  double Psi_bar, Psi_prime_bar;
  afterglow_kernel_Psi(r_bar, &Psi_bar);
  afterglow_kernel_Psi_prime(r_bar, &Psi_prime_bar);

  double delta_r   = r_bar * (delta_c - delta_X);
  double delta_Psi = Psi_prime_bar * delta_r;

  /* Rest-frame sound speed: physical propagation speed in the hidden
     YM sector. For w < 0 fluids, using w directly as the sound speed
     in the pressure-gradient term leads to Jeans-type instability;
     cs2_X (typically 1, causal) keeps perturbations well-posed.
     Analogous to cs2_fld in CLASS's stock fluid module. */
  double cs2 = ap->cs2_X;

  /* (P1) density-contrast equation.
     Standard fluid form with cs2 stabilization:
       delta' = -(1+w)(theta+h'/2) - 3H(cs2-w)delta - 9(1+w)(cs2-ca2)H^2 theta/k^2
     For constant w, ca2 = w so cs2-ca2 = cs2-w.
     Plus MIS non-adiabatic departure: w*(sigma_hat - delta_X). */
  *d_delta_X =
      -one_plus_w * (theta_X + 0.5 * h_prime_syn)
      - 3.0 * H_conf * (cs2 - w_X) * delta_X
      - 9.0 * one_plus_w * (cs2 - w_X) * H_conf * H_conf * theta_X / (k * k)
      - 3.0 * H_conf * w_X * (sigma_hat - delta_X)
      + H_conf * (beta / c_D) * delta_Psi;

  /* (P2) velocity-divergence equation (Euler).
     Pressure gradient uses cs2 (rest-frame sound speed, not w_X). */
  *d_theta_X =
      -H_conf * (1.0 - 3.0 * cs2) * theta_X
      + (k * k) * (cs2 / one_plus_w) * sigma_hat;

  /* (P3) memory equation — the Phase 2 novelty. */
  *d_sigma_hat =
      (H_conf / c_D) * beta * r_bar * Psi_prime_bar * (delta_c - delta_X)
      - (theta_X + 0.5 * h_prime_syn) / (3.0 * c_D)
          * (1.0 - beta * Psi_bar);

  return _SUCCESS_;
}

/* ------------------------------------------------------------------ */
/*  Pressure-perturbation closure: delta p_X = w_X * rho_X_bg * sigma  */
/*                                                                     */
/*  This is the MIS bulk-stress closure. Because delta p_X is tied    */
/*  to sigma_hat (delta Sigma) and NOT to delta_X (delta rho_X), the  */
/*  effective sound speed delta p_X / delta rho_X = w_X * sigma/delta */
/*  is NOT a constant — it becomes w_X only when sigma_hat = delta_X  */
/*  (the adiabatic sub-manifold).                                      */
/* ------------------------------------------------------------------ */
int afterglow_pert_pressure(struct afterglow_params * ap,
                            double rho_X_bg, double delta_X,
                            double sigma_hat,
                            double * delta_p_X) {
  double w_X = -1.0 + 1.0 / (3.0 * ap->c_D);
  double cs2 = ap->cs2_X;
  /* delta_p = cs2 * rho * delta_X + w * rho * (sigma_hat - delta_X)
     = (cs2 - w) * rho * delta_X + w * rho * sigma_hat
     This matches the cs2-stabilized continuity equation form. */
  *delta_p_X = cs2 * rho_X_bg * delta_X + w_X * rho_X_bg * (sigma_hat - delta_X);
  return _SUCCESS_;
}

/* ------------------------------------------------------------------ */
/*  Adiabatic initial conditions                                       */
/*                                                                     */
/*  Deep in the radiation era all fluids share the adiabatic mode:    */
/*                                                                     */
/*    delta_i / (1 + w_i) = delta_gamma / (1 + w_gamma)                */
/*                                                                     */
/*  with w_gamma = 1/3, so 1 + w_gamma = 4/3. For the afterglow:      */
/*                                                                     */
/*    delta_X_ini = (3/4) * (1 + w_X) * delta_gamma                   */
/*                = (3/4) * (1/(3 c_D)) * delta_gamma                 */
/*                = delta_gamma / (4 c_D)                              */
/*                                                                     */
/*  The MIS memory is initialized to the barotropic sub-manifold:     */
/*    sigma_hat_ini = delta_X_ini                                      */
/*  so pressure perturbation delta p_X = w_X * rho_X_bg * delta_X     */
/*  matches the adiabatic sound speed AT IC TIME. Later departures   */
/*  from sigma_hat = delta_X are the direct signature of MIS memory. */
/*                                                                     */
/*  theta_X_ini uses the standard adiabatic relation                  */
/*    theta_i = theta_gamma   (all species co-move initially).        */
/* ------------------------------------------------------------------ */
int afterglow_pert_initial_conditions(struct afterglow_params * ap,
                                      double delta_gamma_IC,
                                      double theta_gamma_IC,
                                      double * delta_X_ini,
                                      double * theta_X_ini,
                                      double * sigma_hat_ini) {
  double w_X = -1.0 + 1.0 / (3.0 * ap->c_D);
  double one_plus_w = 1.0 + w_X;   /* = 1/(3 c_D) */

  *delta_X_ini   = 0.75 * one_plus_w * delta_gamma_IC;  /* adiabatic */
  *theta_X_ini   = theta_gamma_IC;                       /* co-moving */
  *sigma_hat_ini = *delta_X_ini;                          /* barotropic sub-manifold */

  return _SUCCESS_;
}
