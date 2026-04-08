/** @file afterglow.c
 *
 *  Hidden SU(2)_D Yang-Mills afterglow dark-energy module — implementation.
 *
 *  PHASE 1: background sector.
 *  -----------------------------------------------------------------------
 *    Self-contained background integrator for the three new dark-sector
 *    degrees of freedom (rho_X, Sigma, rho_dr,D). This file does NOT yet
 *    modify CLASS's background.c; it is unit-testable in isolation via
 *    test/test_afterglow_bg.py. Phase 1b will wire it into
 *    background_solve_tau().
 *
 *  Reference (all equation numbers verified against the PDF):
 *    Martin & Koh (April 2026), "Dark Energy as the Thermodynamic Afterglow
 *    of a Hidden Gauge-Theory Transition."
 *
 *  Key equations implemented here:
 *
 *    §4.1  Sourced tail on the arrested branch
 *       (30)  rho_X = 3 c_D Sigma,    p_X = -(3 c_D - 1) Sigma
 *       (31)  w_X   = p_X / rho_X  = -1 + 1 / (3 c_D)       [intrinsic EoS]
 *
 *    §4.2  Signed loading kernel
 *       (32)  r     = rho_c / rho_X
 *       (35)  Psi(r) = 4 r (r - 1) / (1 + r)^3
 *       (36)  Psi(r<1) < 0,   Psi(1) = 0,   Psi(r>1) > 0
 *              (extrema at r = 2 ± sqrt(3),  |Psi|_max = 2 sqrt(3)/9)
 *
 *    §4.3  Open-system dynamics (the three RHS equations below)
 *       (37)  u^mu d_mu Sigma = -(Theta/(3 c_D)) [1 - beta Psi(r)] Sigma
 *       (38)  Sigma_dot       = -(H / c_D) [1 - beta Psi(r)] Sigma
 *       (40)  rho_X_dot + 3H (rho_X + p_X) = -Q
 *       (42)  Q               = -(beta / c_D) H Psi(r) rho_X
 *       (43)  rho_X_dot       = -(H / c_D) [1 - beta Psi(r)] rho_X
 *       (44)  rho_c_dot + 3H rho_c = -(beta / c_D) H Psi(r) rho_X
 *
 *    §5.1  Late-time scaling (verified analytically by test 1)
 *       (55)  rho_X ∝ a^(-1/c_D),   w_eff -> w_X = -1 + 1/(3 c_D)
 *
 *  Analytic limits verified by test/test_afterglow_bg.py:
 *
 *    Test 1  beta = 0                  rho_X(N) = Omega0_X * exp(-N/c_D)  [Eq. 43 -> Eq. 55]
 *    Test 2  c_D -> inf, beta = 0      rho_X, Sigma ~ const               [LCDM limit of 38/43]
 *    Test 3  Psi(r) extremum           |Psi|_max = 2 sqrt(3)/9            [Eq. 35]
 *    Test 4  Sign of Psi               Psi(<1)<0, Psi(1)=0, Psi(>1)>0    [Eq. 36]
 *    Test 5  Exchange law Q            sign(Q) = -sign(Psi) for beta>0    [Eq. 42]
 *    Test 6  Smoothness at beta=2.5    finiteness, positivity, monotone   [Eqs. 37 + 42]
 *    Test 7  RHS point-match to Eq.43  max rel err <= 1e-14                [Eq. 43 direct]
 */

#include <math.h>
#include "afterglow/afterglow.h"

/* ------------------------------------------------------------------ */
/*  Signed loading kernel Psi(r)                 [Eq. 35 VERIFIED]    */
/*                                                                     */
/*       Psi(r) = 4 r (r - 1) / (1 + r)^3                             */
/*                                                                     */
/*  Sign convention (Eq. 36):                                         */
/*     r < 1  (matter sub-dominant)  =>  Psi < 0   (stiffening)       */
/*     r = 1  (equality crossing)    =>  Psi = 0                      */
/*     r > 1  (matter dominant)      =>  Psi > 0   (loosening)        */
/*                                                                     */
/*  Extrema at r = 2 ± sqrt(3),  |Psi|_max = 2 sqrt(3)/9 ≈ 0.38490.   */
/* ------------------------------------------------------------------ */
int afterglow_kernel_Psi(double r, double * Psi) {
  double one_plus_r = 1.0 + r;
  double denom = one_plus_r * one_plus_r * one_plus_r;
  *Psi = 4.0 * r * (r - 1.0) / denom;
  return _SUCCESS_;
}

/* ------------------------------------------------------------------ */
/*  Derived exchange law Q                       [Eq. 42 VERIFIED]    */
/*                                                                     */
/*       Q = -(beta / c_D) * H * Psi(r) * rho_X                       */
/*                                                                     */
/*  Follows from combining (38), (40) and the homogeneous relation    */
/*  (29): ρ_conf = V_{D,eq} + 3 c_D Σ.                                 */
/* ------------------------------------------------------------------ */
int afterglow_exchange_Q(double H, double r, double rho_X,
                          struct afterglow_params * ap, double * Q) {
  double Psi;
  afterglow_kernel_Psi(r, &Psi);
  *Q = -(ap->beta_aft / ap->c_D) * H * Psi * rho_X;
  return _SUCCESS_;
}

/* ------------------------------------------------------------------ */
/*  Background right-hand side in N = ln(a)                           */
/*                                                                     */
/*    d rho_X / dN    = -(1/c_D) [1 - beta Psi(r)] rho_X  [Eq. 43]    */
/*    d Sigma / dN    = -(1/c_D) [1 - beta Psi(r)] Sigma  [Eq. 38]    */
/*    d rho_dr,D / dN = -4 rho_dr,D                    (pre-conf rad) */
/*                                                                     */
/*  Equivalently, writing the EoS w_X = -1 + 1/(3 c_D) of Eq. 31,     */
/*                                                                     */
/*    d rho_X / dN    = -3 (1 + w_X) rho_X + (beta / c_D) Psi rho_X   */
/*                                                                     */
/*  The second form is the one that makes contact with the standard   */
/*  CLASS continuity equation  ρ' + 3 H (1+w) ρ = -Q.                 */
/*                                                                     */
/*  [VERIFIED against the afterglow paper, Section 4.3, Eqs. 38–43.]  */
/*                                                                     */
/*  Cold matter in Phase 1 is the undisturbed background              */
/*  Omega0_c * exp(-3 N).  CDM back-reaction from Q (Eq. 44) is       */
/*  activated only in Phase 1b, when CLASS's own background solver    */
/*  is the driver and rho_c is read from the state vector.            */
/* ------------------------------------------------------------------ */
int afterglow_background_rhs(double N, const double * y, double * dy,
                              struct afterglow_params * ap, double Omega0_c) {
  double rho_X = y[AFTERGLOW_IDX_RHO_X];
  double Sigma = y[AFTERGLOW_IDX_SIGMA];
  double rho_dr = y[AFTERGLOW_IDX_RHO_DR_D];

  /* ρ_c(N) — Phase 1: no CDM back-reaction */
  double rho_c = Omega0_c * exp(-3.0 * N);

  /* r = ρ_c / ρ_X   (Eq. 32), guarded */
  double r = (rho_X > 1e-300) ? (rho_c / rho_X) : 0.0;

  double Psi;
  afterglow_kernel_Psi(r, &Psi);

  /* d ρ_X / dN = -(1/c_D)[1 - β Ψ(r)] ρ_X                  (Eq. 43) */
  dy[AFTERGLOW_IDX_RHO_X] =
      -(1.0 / ap->c_D) * (1.0 - ap->beta_aft * Psi) * rho_X;

  /* d Σ  / dN = -(1/c_D)[1 - β Ψ(r)] Σ                     (Eq. 38) */
  dy[AFTERGLOW_IDX_SIGMA] =
      -(1.0 / ap->c_D) * (1.0 - ap->beta_aft * Psi) * Sigma;

  /* d ρ_dr,D / dN = -4 ρ_dr,D   (free-streaming radiation before
     confinement; post-confinement transfer Q_c at z = z_conf will
     be handled in Phase 1b when it is known which side of z_conf
     the stepper is on)                                              */
  dy[AFTERGLOW_IDX_RHO_DR_D] = -4.0 * rho_dr;

  return _SUCCESS_;
}

/* ------------------------------------------------------------------ */
/*  4th-order Runge-Kutta integrator                                   */
/*                                                                     */
/*  The grid N_grid[] is monotonically decreasing (N_grid[0] = 0      */
/*  today, later entries are more negative). Each step of size         */
/*  dN_i = N_grid[i] - N_grid[i-1] (negative) is integrated via RK4.  */
/* ------------------------------------------------------------------ */
int afterglow_background_evolve(struct afterglow_params * ap,
                                 double Omega0_X, double Omega0_c,
                                 int N_steps, const double * N_grid,
                                 double * rho_X_out, double * Sigma_out,
                                 double * rho_dr_out) {
  double y[AFTERGLOW_NBG];
  double k1[AFTERGLOW_NBG], k2[AFTERGLOW_NBG];
  double k3[AFTERGLOW_NBG], k4[AFTERGLOW_NBG];
  double ytmp[AFTERGLOW_NBG];
  int i, j;

  if (N_steps < 2) return _FAILURE_;

  /* initial conditions at N = 0: ρ_X(0) = Ω₀_X,  Σ(0) = Ω₀_X/(3 c_D) by Eq. 30 */
  y[AFTERGLOW_IDX_RHO_X] = Omega0_X;
  y[AFTERGLOW_IDX_SIGMA] = (ap->Sigma_today >= 0.0)
                           ? ap->Sigma_today
                           : (Omega0_X / (3.0 * ap->c_D));
  y[AFTERGLOW_IDX_RHO_DR_D] = ap->Omega0_dr_dark;

  rho_X_out[0]  = y[AFTERGLOW_IDX_RHO_X];
  Sigma_out[0]  = y[AFTERGLOW_IDX_SIGMA];
  rho_dr_out[0] = y[AFTERGLOW_IDX_RHO_DR_D];

  for (i = 1; i < N_steps; i++) {
    double N  = N_grid[i - 1];
    double dN = N_grid[i] - N_grid[i - 1];
    double h2 = 0.5 * dN;

    afterglow_background_rhs(N, y, k1, ap, Omega0_c);
    for (j = 0; j < AFTERGLOW_NBG; j++) ytmp[j] = y[j] + h2 * k1[j];
    afterglow_background_rhs(N + h2, ytmp, k2, ap, Omega0_c);
    for (j = 0; j < AFTERGLOW_NBG; j++) ytmp[j] = y[j] + h2 * k2[j];
    afterglow_background_rhs(N + h2, ytmp, k3, ap, Omega0_c);
    for (j = 0; j < AFTERGLOW_NBG; j++) ytmp[j] = y[j] + dN * k3[j];
    afterglow_background_rhs(N + dN, ytmp, k4, ap, Omega0_c);

    for (j = 0; j < AFTERGLOW_NBG; j++) {
      y[j] += (dN / 6.0) * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
    }

    rho_X_out[i]  = y[AFTERGLOW_IDX_RHO_X];
    Sigma_out[i]  = y[AFTERGLOW_IDX_SIGMA];
    rho_dr_out[i] = y[AFTERGLOW_IDX_RHO_DR_D];
  }

  return _SUCCESS_;
}

/* ------------------------------------------------------------------ */
/*  Lifecycle                                                          */
/* ------------------------------------------------------------------ */
int afterglow_init(struct afterglow_params * ap) {
  (void) ap;
  return _SUCCESS_;
}

int afterglow_free(struct afterglow_params * ap) {
  (void) ap;
  return _SUCCESS_;
}
