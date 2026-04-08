/** @file afterglow.h
 *
 *  Hidden SU(2)_D Yang-Mills afterglow dark-energy module for CLASS_SYMT.
 *
 *  Implements the causal imperfect fluid described in Martin & Koh (2026),
 *  "Dark Energy as the Thermodynamic Afterglow of a Hidden Gauge-Theory
 *  Transition." The fluid carries a Müller-Israel-Stewart memory variable
 *  Sigma(tau) and couples to the cold matter sector via the signed kernel
 *
 *      Psi(r) = 4 r (r - 1) / (1 + r)^3,   r = rho_c / rho_X
 *
 *  The loading-unloading law (Eq. 37) and derived exchange law (Eq. 42) are:
 *
 *      u^mu d_mu Sigma = -(Theta / (3 c_D)) * [1 - beta * Psi(r)] * Sigma
 *      Q               = -(beta / c_D) * H * Psi(r) * rho_X
 *
 *  PHASE 0 SCAFFOLD: prototypes and struct only. Integrator stubs live in
 *  source/afterglow/afterglow.c and currently return _SUCCESS_ without
 *  modifying CLASS state (LCDM limit).
 */

#ifndef __AFTERGLOW__
#define __AFTERGLOW__

#include "common.h"

/* ─── User parameters (parsed from .ini) ─────────────────────────────── */
struct afterglow_params {
  double Omega0_dr_dark;    /**< pre-confinement dark radiation today     */
  double z_conf;            /**< confinement redshift                     */
  double c_D;               /**< Deborah number (relaxation / Hubble)     */
  double beta_aft;          /**< loading coupling beta (Eq. 37, Eq. 42)   */
  double Sigma_today;       /**< IC for memory variable; <0 means auto    */
  double cs2_X;             /**< rest-frame sound speed^2 (causal => 1)   */
  short afterglow_on;       /**< 0 => LCDM limit, 1 => full model         */
};

/* ─── Background indices (filled by afterglow_init) ──────────────────── */
struct afterglow_bg_indices {
  int index_bg_rho_X;
  int index_bg_Sigma;
  int index_bg_rho_dr_dark;
  int index_bg_p_X;
  int bg_size;
};

/* ─── Perturbation indices (Phase 2) ─────────────────────────────────── */
struct afterglow_pt_indices {
  int index_pt_delta_X;
  int index_pt_theta_X;
  int index_pt_delta_Sigma;   /* CRITICAL: independent DOF, do not fold   */
  int pt_size;
};

/* ─── API (scaffold — implementations return _SUCCESS_) ──────────────── */
int afterglow_kernel_Psi(double r, double * Psi);
int afterglow_background_derivs(double z, double * y, double * dy,
                                 struct afterglow_params * ap);
int afterglow_exchange_Q(double H, double r, double rho_X,
                          struct afterglow_params * ap, double * Q);
int afterglow_init(struct afterglow_params * ap);
int afterglow_free(struct afterglow_params * ap);

#endif /* __AFTERGLOW__ */
