/** @file afterglow.h
 *
 *  Hidden SU(2)_D Yang-Mills afterglow dark-energy module for CLASS_SYMT.
 *
 *  Implements the causal imperfect fluid described in Martin & Koh (April 2026),
 *  "Dark Energy as the Thermodynamic Afterglow of a Hidden Gauge-Theory
 *  Transition." The fluid carries a Müller-Israel-Stewart memory variable
 *  Sigma and couples to the cold matter sector via the signed kernel
 *
 *      (35)  Psi(r) = 4 r (r - 1) / (1 + r)^3,   r = rho_c / rho_X   [§4.2]
 *
 *  Paper equations used below (ALL VERIFIED against the PDF):
 *
 *      (31)  w_X = p_X / rho_X = -1 + 1/(3 c_D)                       [§4.1]
 *      (37)  u^mu d_mu Sigma = -(Theta/(3 c_D)) [1 - beta Psi(r)] Sigma [§4.3]
 *      (38)  Sigma_dot       = -(H / c_D) [1 - beta Psi(r)] Sigma     [§4.3]
 *      (42)  Q               = -(beta / c_D) H Psi(r) rho_X           [§4.3]
 *      (43)  rho_X_dot       = -(H / c_D) [1 - beta Psi(r)] rho_X     [§4.3]
 *      (55)  rho_X ∝ a^(-1/c_D),  w_eff -> w_X                         [§5.1]
 *
 *  PHASE 1: background sector complete.
 *    - afterglow_kernel_Psi()            : signed kernel Psi(r)
 *    - afterglow_exchange_Q()            : derived exchange law Q
 *    - afterglow_background_rhs()        : d/dN of {rho_X, Sigma, rho_dr_D}
 *    - afterglow_background_evolve()     : 4th-order Runge-Kutta integrator
 *                                          over a user-supplied N = ln(a) grid
 *
 *  Units: energy densities are carried in units of rho_crit_0 so that
 *  rho_X(N=0) = Omega0_X at a = 1. The caller multiplies by rho_crit_0
 *  before handing the values to the CLASS background state vector.
 */

#ifndef __AFTERGLOW__
#define __AFTERGLOW__

#include "common.h"

/* ─── Derivative vector layout ───────────────────────────────────────── */
#define AFTERGLOW_NBG 3             /**< number of background DOFs        */
#define AFTERGLOW_IDX_RHO_X     0   /**< index for rho_X / rho_crit_0     */
#define AFTERGLOW_IDX_SIGMA     1   /**< index for Sigma / rho_crit_0     */
#define AFTERGLOW_IDX_RHO_DR_D  2   /**< index for rho_dr,D / rho_crit_0  */

/* ─── User parameters (parsed from .ini) ─────────────────────────────── */
struct afterglow_params {
  double Omega0_dr_dark;    /**< pre-confinement dark radiation today     */
  double z_conf;            /**< confinement redshift                     */
  double c_D;               /**< Deborah number (relaxation / Hubble)     */
  double beta_aft;          /**< loading coupling beta (Eq. 37, Eq. 42)   */
  double Sigma_today;       /**< IC for memory variable; <0 means auto    */
  double cs2_X;             /**< rest-frame sound speed^2 (causal => 1)   */
  double bp_regulator;      /**< B1 regulator: smooth max(0, 1-beta*Psi)
                                 with width delta. 0 => disabled (default;
                                 preserves Phase 2 regression). >0 enables
                                 smooth clamp in afterglow_pert_rhs P3 to
                                 prevent unstable mode when beta*Psi > 1.   */
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

/* ─── Core physics ───────────────────────────────────────────────────── */

/**
 *  Signed loading kernel  Psi(r) = 4 r (r - 1) / (1 + r)^3.  Pure function.
 */
int afterglow_kernel_Psi(double r, double * Psi);

/**
 *  Derived exchange law  Q = -(beta / c_D) * H * Psi(r) * rho_X  (Eq. 42).
 */
int afterglow_exchange_Q(double H, double r, double rho_X,
                          struct afterglow_params * ap, double * Q);

/**
 *  Right-hand side d/dN of the background vector y = {rho_X, Sigma, rho_dr_D}.
 *
 *  N = ln(a), with a = 1 today. Equations integrated:
 *
 *      d rho_X / dN     = -(1/c_D) [1 - beta * Psi(r)] * rho_X           (Eq. 43)
 *      d Sigma / dN     = -(1/c_D) [1 - beta * Psi(r)] * Sigma            (Eq. 38)
 *      d rho_dr,D / dN  = -4 rho_dr,D                                      (free radiation)
 *
 *  with w_X = -1 + 1/(3 c_D) [Eq. 31] and r = rho_c / rho_X [Eq. 32].
 *  For Phase 1 the cold
 *  matter is  Ω_c0 * exp(-3N)  (no CDM backreaction yet).
 */
int afterglow_background_rhs(double N, const double * y, double * dy,
                              struct afterglow_params * ap, double Omega0_c);

/**
 *  4th-order Runge-Kutta integrator for the afterglow background sector.
 *  Evolves backward in N = ln(a) from N_grid[0] = 0 (today) to N_grid[N_steps-1].
 */
int afterglow_background_evolve(struct afterglow_params * ap,
                                 double Omega0_X, double Omega0_c,
                                 int N_steps, const double * N_grid,
                                 double * rho_X_out, double * Sigma_out,
                                 double * rho_dr_out);

/* ─── Lifecycle ──────────────────────────────────────────────────────── */
int afterglow_init(struct afterglow_params * ap);
int afterglow_free(struct afterglow_params * ap);

#endif /* __AFTERGLOW__ */
