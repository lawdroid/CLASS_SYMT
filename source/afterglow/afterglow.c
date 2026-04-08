/** @file afterglow.c
 *
 *  Hidden SU(2)_D Yang-Mills afterglow dark-energy module — implementation.
 *  Phase 0 scaffold: the entry points exist and return _SUCCESS_ so that the
 *  rest of CLASS runs unmodified (LCDM limit). Phase 1 fills in the background
 *  ODE integrators on branch `feature/afterglow-background`.
 *
 *  Reference: Martin & Koh (2026), "Dark Energy as the Thermodynamic
 *  Afterglow of a Hidden Gauge-Theory Transition."
 */

#include "afterglow/afterglow.h"

/**
 * Signed loading kernel Psi(r) = 4 r (r - 1) / (1 + r)^3.
 *
 * Sign convention:
 *   r < 1  (matter subdominant)  => Psi < 0  (stiffening)
 *   r = 1  (equality crossing)   => Psi = 0
 *   r > 1  (matter dominant)     => Psi > 0  (loosening)
 *
 * Max |Psi| = 2 sqrt(3) / 9 ≈ 0.3849 at r = 2 ± sqrt(3).
 */
int afterglow_kernel_Psi(double r, double * Psi) {
  double one_plus_r = 1.0 + r;
  *Psi = 4.0 * r * (r - 1.0) / (one_plus_r * one_plus_r * one_plus_r);
  return _SUCCESS_;
}

/**
 * Derived exchange law (Eq. 42): Q = -(beta / c_D) * H * Psi(r) * rho_X.
 */
int afterglow_exchange_Q(double H, double r, double rho_X,
                          struct afterglow_params * ap, double * Q) {
  double Psi;
  afterglow_kernel_Psi(r, &Psi);
  *Q = -(ap->beta_aft / ap->c_D) * H * Psi * rho_X;
  return _SUCCESS_;
}

/**
 * Background derivative vector (stub).
 *
 * PHASE 1 TODO: populate dy[index_bg_rho_X], dy[index_bg_Sigma],
 * dy[index_bg_rho_dr_dark] with
 *
 *   d rho_X / dN      = -3 (1 + w_X) rho_X - beta * Psi(r) * rho_X
 *   d Sigma / dN      = -[1 - beta * Psi(r)] * Sigma / c_D
 *   d rho_dr_D / dN   = -4 rho_dr_D - Q_c(z; z_conf)
 *
 * where N = ln a and Q_c is a smooth tanh activating at z_conf.
 */
int afterglow_background_derivs(double z, double * y, double * dy,
                                 struct afterglow_params * ap) {
  /* scaffold: no-op, see Phase 1 TODO above */
  (void) z; (void) y; (void) dy; (void) ap;
  return _SUCCESS_;
}

int afterglow_init(struct afterglow_params * ap) {
  /* Phase 0: no allocations yet */
  (void) ap;
  return _SUCCESS_;
}

int afterglow_free(struct afterglow_params * ap) {
  (void) ap;
  return _SUCCESS_;
}
