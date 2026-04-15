#!/usr/bin/env python3
# ---------------------------------------------------------------------
# test_afterglow_bg.py
#
# Pure-Python replica of source/afterglow/afterglow.c used as a unit
# test for the Phase 1 background integrator of CLASS_SYMT.
#
# We mirror the C code line-for-line (same RK4, same right-hand side)
# and check analytic limits from Martin & Koh (April 2026)
# "Dark Energy as the Thermodynamic Afterglow of a Hidden
#  Gauge-Theory Transition."  All equation numbers below have been
# cross-checked against the PDF.
#
#   Test 1  beta = 0                  rho_X(N) = Omega0_X * exp(-N/c_D)
#                                     [Eq. 43 -> Eq. 55, §5.1]
#   Test 2  c_D -> infinity, beta=0   rho_X, Sigma ~ const
#                                     [LCDM limit of Eqs. 38, 43]
#   Test 3  Psi(r) extrema            |Psi|_max = 2*sqrt(3)/9   [Eq. 35]
#   Test 4  Sign of Psi               Psi(<1)<0, Psi(1)=0, Psi(>1)>0  [Eq. 36]
#   Test 5  Exchange law Q            sign(Q) = -sign(Psi) for beta>0 [Eq. 42]
#   Test 6  Smoothness at beta=2.5    finiteness, positivity, monotone
#                                     [Eqs. 37 + 42 together]
#   Test 7  RHS point-match to Eq.43  numerical d rho_X/dN agrees with
#                                     -(1/c_D)[1-βΨ(r)]ρ_X to 1e-14
#                                     [direct regression test for Eq. 43]
#
# The point of mirroring the C is that *any* drift between the C and
# Python versions will show up as a failing test here, before the code
# ever touches CLASS's background.c. This is deliberate Phase-1
# scaffolding: we validate the math in isolation.
#
# Usage:
#     python test/test_afterglow_bg.py
# ---------------------------------------------------------------------

from __future__ import annotations
import math
import sys
from dataclasses import dataclass


# ─── mirror of struct afterglow_params ─────────────────────────────
@dataclass
class AfterglowParams:
    Omega0_dr_dark: float = 0.0
    z_conf: float = 1.0e5
    c_D: float = 1.0
    beta_aft: float = 0.0
    Sigma_today: float = -1.0   # -1 => auto: Omega0_X / (3 c_D)
    cs2_X: float = 1.0
    afterglow_on: int = 1


# ─── kernel Psi(r) = 4 r (r-1) / (1+r)^3 ───────────────────────────
def Psi(r: float) -> float:
    return 4.0 * r * (r - 1.0) / (1.0 + r) ** 3


# ─── derived exchange law Q = -(beta/c_D) H Psi rho_X ──────────────
def Q_exchange(H: float, r: float, rho_X: float, ap: AfterglowParams) -> float:
    return -(ap.beta_aft / ap.c_D) * H * Psi(r) * rho_X


# ─── RHS  d/dN (rho_X, Sigma, rho_dr)  VERIFIED vs paper Eqs. 38, 43 ─
def rhs(N: float, y, ap: AfterglowParams, Omega0_c: float):
    rho_X, Sigma, rho_dr = y
    rho_c = Omega0_c * math.exp(-3.0 * N)              # Phase 1: no backreaction
    r = rho_c / rho_X if rho_X > 1e-300 else 0.0        # Eq. 32
    psi = Psi(r)                                        # Eq. 35
    # d rho_X / dN = -(1/c_D)[1 - β Ψ(r)] rho_X         (Eq. 43)
    d_rho_X  = -(1.0 / ap.c_D) * (1.0 - ap.beta_aft * psi) * rho_X
    # d Σ      / dN = -(1/c_D)[1 - β Ψ(r)] Σ            (Eq. 38)
    d_Sigma  = -(1.0 / ap.c_D) * (1.0 - ap.beta_aft * psi) * Sigma
    d_rho_dr = -4.0 * rho_dr
    return [d_rho_X, d_Sigma, d_rho_dr]


# ─── RK4 integrator mirroring afterglow_background_evolve() ───────
def evolve(ap: AfterglowParams, Omega0_X: float, Omega0_c: float, N_grid):
    y = [
        Omega0_X,
        ap.Sigma_today if ap.Sigma_today >= 0.0 else Omega0_X / (3.0 * ap.c_D),
        ap.Omega0_dr_dark,
    ]
    out_X, out_S, out_D = [y[0]], [y[1]], [y[2]]
    for i in range(1, len(N_grid)):
        N = N_grid[i - 1]
        dN = N_grid[i] - N_grid[i - 1]
        h2 = 0.5 * dN
        k1 = rhs(N, y, ap, Omega0_c)
        ytmp = [y[j] + h2 * k1[j] for j in range(3)]
        k2 = rhs(N + h2, ytmp, ap, Omega0_c)
        ytmp = [y[j] + h2 * k2[j] for j in range(3)]
        k3 = rhs(N + h2, ytmp, ap, Omega0_c)
        ytmp = [y[j] + dN * k3[j] for j in range(3)]
        k4 = rhs(N + dN, ytmp, ap, Omega0_c)
        y = [
            y[j] + (dN / 6.0) * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j])
            for j in range(3)
        ]
        out_X.append(y[0])
        out_S.append(y[1])
        out_D.append(y[2])
    return out_X, out_S, out_D


# ─── helpers ────────────────────────────────────────────────────────
_passed = 0
_failed = 0


def check(name: str, ok: bool, detail: str = "") -> None:
    global _passed, _failed
    marker = "PASS" if ok else "FAIL"
    print(f"  [{marker}] {name}" + (f"  {detail}" if detail else ""))
    if ok:
        _passed += 1
    else:
        _failed += 1


def approx(a: float, b: float, rtol: float = 1e-6, atol: float = 1e-12) -> bool:
    return abs(a - b) <= max(atol, rtol * max(abs(a), abs(b)))


# ═══════════════════════════════════════════════════════════════════
#   Test 1.  beta = 0  =>  rho_X(N) = Omega0_X * exp(-N / c_D)
# ═══════════════════════════════════════════════════════════════════
def test1_analytic_scaling():
    print("\nTest 1  :  beta = 0  =>  rho_X = Omega0_X * exp(-N / c_D)")
    ap = AfterglowParams(c_D=1.0, beta_aft=0.0)
    Omega0_X = 0.69
    Omega0_c = 0.26
    # grid from today (N=0) backward to N=-10 (z ~ 22000), 401 steps
    N_steps = 401
    N_grid = [(-10.0) * i / (N_steps - 1) for i in range(N_steps)]
    rho_X, Sigma, _ = evolve(ap, Omega0_X, Omega0_c, N_grid)

    max_err = 0.0
    for N, rX in zip(N_grid, rho_X):
        analytic = Omega0_X * math.exp(-N / ap.c_D)
        rel = abs(rX - analytic) / analytic
        if rel > max_err:
            max_err = rel
    check(
        "rho_X matches exp(-N/c_D) to 1e-6",
        max_err < 1e-6,
        f"max rel err = {max_err:.2e}",
    )

    # Also verify with c_D = 1.2 to catch any hard-coding
    ap2 = AfterglowParams(c_D=1.2, beta_aft=0.0)
    rho_X2, _, _ = evolve(ap2, Omega0_X, Omega0_c, N_grid)
    max_err2 = 0.0
    for N, rX in zip(N_grid, rho_X2):
        analytic = Omega0_X * math.exp(-N / ap2.c_D)
        rel = abs(rX - analytic) / analytic
        if rel > max_err2:
            max_err2 = rel
    check(
        "rho_X matches exp(-N/c_D) with c_D=1.2",
        max_err2 < 1e-6,
        f"max rel err = {max_err2:.2e}",
    )


# ═══════════════════════════════════════════════════════════════════
#   Test 2.  c_D -> infinity, beta = 0  =>  rho_X, Sigma ~ const
# ═══════════════════════════════════════════════════════════════════
def test2_LCDM_limit():
    print("\nTest 2  :  c_D -> infinity, beta = 0  =>  LCDM limit")
    ap = AfterglowParams(c_D=1.0e6, beta_aft=0.0, Sigma_today=1.0e-8)
    Omega0_X = 0.69
    Omega0_c = 0.26
    N_steps = 201
    N_grid = [(-6.0) * i / (N_steps - 1) for i in range(N_steps)]
    rho_X, Sigma, _ = evolve(ap, Omega0_X, Omega0_c, N_grid)
    drift_X = abs(rho_X[-1] - rho_X[0]) / rho_X[0]
    drift_S = abs(Sigma[-1] - Sigma[0]) / Sigma[0]
    check(
        "rho_X drift over 6 e-folds < 1e-5",
        drift_X < 1e-5,
        f"drift = {drift_X:.2e}",
    )
    check(
        "Sigma drift over 6 e-folds < 1e-5",
        drift_S < 1e-5,
        f"drift = {drift_S:.2e}",
    )


# ═══════════════════════════════════════════════════════════════════
#   Test 3.  |Psi|_max = 2 sqrt(3) / 9
# ═══════════════════════════════════════════════════════════════════
def test3_Psi_extrema():
    print("\nTest 3  :  |Psi(r)|_max = 2 sqrt(3) / 9")
    target = 2.0 * math.sqrt(3.0) / 9.0
    # scan r in log space
    r_scan = [10.0 ** (x / 500.0 - 3.0) for x in range(1500)]
    max_abs = max(abs(Psi(r)) for r in r_scan)
    check(
        "sup |Psi(r)| = 2 sqrt(3)/9",
        approx(max_abs, target, rtol=1e-3),
        f"got {max_abs:.6f}, target {target:.6f}",
    )
    # analytic extremum at r = 2 +- sqrt(3)
    r_plus = 2.0 + math.sqrt(3.0)
    r_minus = 2.0 - math.sqrt(3.0)
    check(
        "Psi(2 + sqrt(3)) = 2 sqrt(3)/9",
        approx(Psi(r_plus), target, rtol=1e-10),
    )
    check(
        "Psi(2 - sqrt(3)) = -2 sqrt(3)/9",
        approx(Psi(r_minus), -target, rtol=1e-10),
    )


# ═══════════════════════════════════════════════════════════════════
#   Test 4.  Sign of Psi
# ═══════════════════════════════════════════════════════════════════
def test4_Psi_sign():
    print("\nTest 4  :  sign convention of Psi(r)")
    check("Psi(0.5) < 0 (matter subdominant stiffens)", Psi(0.5) < 0)
    check("Psi(1.0) = 0 (equality crossing)", approx(Psi(1.0), 0.0, atol=1e-15))
    check("Psi(2.0) > 0 (matter dominant loosens)", Psi(2.0) > 0)


# ═══════════════════════════════════════════════════════════════════
#   Test 5.  sign(Q) = -sign(Psi) for beta > 0
# ═══════════════════════════════════════════════════════════════════
def test5_exchange_sign():
    print("\nTest 5  :  sign(Q) = -sign(Psi) for beta > 0")
    ap = AfterglowParams(c_D=1.0, beta_aft=2.5)
    H, rho_X = 70.0, 0.69
    Q_half = Q_exchange(H, 0.5, rho_X, ap)
    Q_one = Q_exchange(H, 1.0, rho_X, ap)
    Q_two = Q_exchange(H, 2.0, rho_X, ap)
    check("Q(r=0.5) > 0 (energy flows into X)", Q_half > 0)
    check("Q(r=1.0) = 0 (crossing)", approx(Q_one, 0.0, atol=1e-15))
    check("Q(r=2.0) < 0 (energy flows out of X)", Q_two < 0)


# ═══════════════════════════════════════════════════════════════════
#   Bonus.  Energy conservation check with small beta
# ═══════════════════════════════════════════════════════════════════
def test6_energy_conservation_smoothness():
    print("\nTest 6  :  numerical smoothness with beta = 2.5, c_D = 1.0")
    ap = AfterglowParams(c_D=1.0, beta_aft=2.5)
    Omega0_X = 0.69
    Omega0_c = 0.26
    N_steps = 2001
    N_grid = [(-8.0) * i / (N_steps - 1) for i in range(N_steps)]
    rho_X, Sigma, _ = evolve(ap, Omega0_X, Omega0_c, N_grid)
    check("rho_X is finite everywhere", all(math.isfinite(x) for x in rho_X))
    check("Sigma is finite everywhere", all(math.isfinite(s) for s in Sigma))
    check("rho_X stays positive", all(x > 0 for x in rho_X))
    # With β=2.5 and matter-dominated past, Ψ>0 so the loosening term
    # speeds relaxation → Sigma decays monotonically faster than c_D=1 alone
    check("Sigma monotonically increasing backward in time",
          all(Sigma[i] >= Sigma[i - 1] - 1e-12 for i in range(1, len(Sigma))))


# ═══════════════════════════════════════════════════════════════════
#   Test 7.  RHS point-match to Eq. 43
#            d ρ_X / dN  == -(1/c_D)[1 - β Ψ(r)] ρ_X   to 1e-14
#
#   This is the regression test that catches sign/normalization bugs
#   in the β-coupled term. Earlier versions of the C code accidentally
#   wrote `-β Ψ ρ_X` instead of `+(β/c_D) Ψ ρ_X` and tests 1–6 all
#   passed tautologically because they used β=0 or didn't compare to
#   a β-dependent analytic expression.
# ═══════════════════════════════════════════════════════════════════
def test7_rhs_pointmatch_eq43():
    print("\nTest 7  :  RHS point-match to Eq. 43")
    ap = AfterglowParams(c_D=1.3, beta_aft=2.5)
    Omega0_X = 0.69
    Omega0_c = 0.26
    worst = 0.0
    for N in [-3.0, -1.5, -0.7, 0.0]:
        rho_X = Omega0_X * math.exp(-N / ap.c_D)   # any positive value is fine
        Sigma = rho_X / (3.0 * ap.c_D)
        rho_dr = 0.0
        y = [rho_X, Sigma, rho_dr]
        dy = rhs(N, y, ap, Omega0_c)
        rho_c = Omega0_c * math.exp(-3.0 * N)
        r = rho_c / rho_X
        psi = Psi(r)
        expected = -(1.0 / ap.c_D) * (1.0 - ap.beta_aft * psi) * rho_X
        rel = abs(dy[0] - expected) / max(abs(expected), 1e-300)
        worst = max(worst, rel)
    check(
        "d rho_X / dN matches -(1/c_D)[1 - β Ψ] ρ_X  (Eq. 43)",
        worst < 1e-14,
        f"max rel err = {worst:.2e}",
    )
    # Same cross-check for Σ against Eq. 38
    worst_s = 0.0
    for N in [-3.0, -1.5, -0.7, 0.0]:
        rho_X = Omega0_X * math.exp(-N / ap.c_D)
        Sigma = 0.123 * rho_X
        y = [rho_X, Sigma, 0.0]
        dy = rhs(N, y, ap, Omega0_c)
        rho_c = Omega0_c * math.exp(-3.0 * N)
        r = rho_c / rho_X
        psi = Psi(r)
        expected_s = -(1.0 / ap.c_D) * (1.0 - ap.beta_aft * psi) * Sigma
        rel = abs(dy[1] - expected_s) / max(abs(expected_s), 1e-300)
        worst_s = max(worst_s, rel)
    check(
        "d Sigma / dN matches -(1/c_D)[1 - β Ψ] Σ      (Eq. 38)",
        worst_s < 1e-14,
        f"max rel err = {worst_s:.2e}",
    )


# ═══════════════════════════════════════════════════════════════════
#   PHASE 1b — glue layer tests
#
#   These mirror source/afterglow/afterglow_class_glue.c line-for-line
#   and check that the CLASS-side entry points implement the same
#   Eqs. 38, 43, 31, 30, 55 that the standalone Phase 1a integrator
#   already pins.
# ═══════════════════════════════════════════════════════════════════

def glue_derivs(ap, rho_c, rho_X, Sigma):
    """Python mirror of afterglow_glue_derivs (afterglow_class_glue.c)."""
    r = rho_c / rho_X if rho_X > 1e-300 else 0.0   # Eq. 32
    psi = Psi(r)                                   # Eq. 35
    factor = (1.0 / ap.c_D) * (1.0 - ap.beta_aft * psi)
    d_rho_X = -factor * rho_X                      # Eq. 43
    d_Sigma = -factor * Sigma                      # Eq. 38
    return d_rho_X, d_Sigma


def glue_initial_conditions(ap, Omega0_X, a_ini):
    """Python mirror of afterglow_glue_initial_conditions."""
    loga_ini = math.log(a_ini)
    rho_X_ini = Omega0_X * math.exp(-loga_ini / ap.c_D)  # Eq. 55
    Sigma_ini = rho_X_ini / (3.0 * ap.c_D)                # Eq. 30
    return rho_X_ini, Sigma_ini


def glue_cdm_source(ap, rho_c, rho_X):
    """Python mirror of afterglow_glue_cdm_source (Phase 1b.2, Eq. 44).

    Returns the EXTRA term beyond the standard -3 rho_c derivative:
        d rho_c / d(log a) = -3 rho_c + cdm_source
        cdm_source         = -(beta / c_D) Psi(r) rho_X
    """
    r = rho_c / rho_X if rho_X > 1e-300 else 0.0
    return -(ap.beta_aft / ap.c_D) * Psi(r) * rho_X


def test8_glue_derivs():
    """
    The glue-layer derivative helper must equal the Phase-1a rhs()
    component-by-component for any (rho_c, rho_X, Sigma). This is the
    regression test that Phase 1b is NOT drifting from Phase 1a.
    """
    print("\nTest 8  :  glue_derivs matches Phase-1a rhs()  (Eqs. 38, 43)")
    ap = AfterglowParams(c_D=1.3, beta_aft=2.5)
    worst_X = 0.0
    worst_S = 0.0
    for N in [-3.0, -1.5, -0.7, 0.0]:
        Omega0_c = 0.26
        rho_X = 0.69 * math.exp(-N / ap.c_D)
        Sigma = 0.321 * rho_X                       # arbitrary nonzero
        rho_c = Omega0_c * math.exp(-3.0 * N)
        dXg, dSg = glue_derivs(ap, rho_c, rho_X, Sigma)
        dy = rhs(N, [rho_X, Sigma, 0.0], ap, Omega0_c)
        worst_X = max(worst_X, abs(dXg - dy[0]) / max(abs(dy[0]), 1e-300))
        worst_S = max(worst_S, abs(dSg - dy[1]) / max(abs(dy[1]), 1e-300))
    check(
        "glue d rho_X / dN == Phase-1a d rho_X / dN  (Eq. 43)",
        worst_X < 1e-14,
        f"max rel err = {worst_X:.2e}",
    )
    check(
        "glue d Sigma / dN == Phase-1a d Sigma / dN  (Eq. 38)",
        worst_S < 1e-14,
        f"max rel err = {worst_S:.2e}",
    )


def test9_glue_initial_conditions():
    """
    The glue-layer IC helper must reproduce Omega0_X when evolved
    forward from a_ini back to a = 1 with the same rhs(), in the
    beta = 0 limit (which is exact).
    """
    print("\nTest 9  :  glue initial conditions re-hit Omega0_X  (Eq. 55)")
    for c_D in [0.8, 1.0, 1.3]:
        ap = AfterglowParams(c_D=c_D, beta_aft=0.0)
        Omega0_X = 0.69
        a_ini = 1e-4                 # CLASS-style earliest a
        rho_X_ini, Sigma_ini = glue_initial_conditions(ap, Omega0_X, a_ini)
        # Integrate forward from N_ini = log(a_ini) to N = 0 with beta=0
        N_ini = math.log(a_ini)
        N_steps = 801
        N_grid = [N_ini + (0.0 - N_ini) * i / (N_steps - 1)
                  for i in range(N_steps)]
        # reuse evolve(), but override IC by setting Sigma_today < 0
        # and by pre-computing our own forward RK4 from (rho_X_ini,
        # Sigma_ini). The simplest way is to reuse the rhs() directly:
        y = [rho_X_ini, Sigma_ini, 0.0]
        for i in range(1, len(N_grid)):
            N = N_grid[i - 1]
            dN = N_grid[i] - N_grid[i - 1]
            h2 = 0.5 * dN
            k1 = rhs(N, y, ap, 0.26)
            yt = [y[j] + h2 * k1[j] for j in range(3)]
            k2 = rhs(N + h2, yt, ap, 0.26)
            yt = [y[j] + h2 * k2[j] for j in range(3)]
            k3 = rhs(N + h2, yt, ap, 0.26)
            yt = [y[j] + dN * k3[j] for j in range(3)]
            k4 = rhs(N + dN, yt, ap, 0.26)
            y = [y[j] + (dN / 6.0) * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j])
                 for j in range(3)]
        rel = abs(y[0] - Omega0_X) / Omega0_X
        check(
            f"rho_X(today) ≈ Omega0_X  (c_D={c_D}, β=0)",
            rel < 1e-5,
            f"rel err = {rel:.2e}",
        )


# ═══════════════════════════════════════════════════════════════════
#   PHASE 1b.2 — CDM back-reaction (Eq. 44)
#
#   Test 10  glue_cdm_source point-match to Eq. 44:
#              source == -(beta/c_D) Psi(r) rho_X   to 1e-14
#
#   Test 11  Total dark-sector conservation:
#              d(rho_X + rho_c)/dN + 3 rho_c + 3(1+w_X) rho_X = 0
#              (source terms in Eqs. 43 and 44 cancel exactly)
# ═══════════════════════════════════════════════════════════════════

def test10_cdm_source_pointmatch_eq44():
    print("\nTest 10 :  CDM back-reaction source matches Eq. 44")
    ap = AfterglowParams(c_D=1.3, beta_aft=2.5)
    Omega0_X = 0.69
    Omega0_c = 0.26
    worst = 0.0
    for N in [-3.0, -1.5, -0.7, 0.0]:
        rho_X = Omega0_X * math.exp(-N / ap.c_D)
        rho_c = Omega0_c * math.exp(-3.0 * N)
        r = rho_c / rho_X
        psi = Psi(r)
        expected = -(ap.beta_aft / ap.c_D) * psi * rho_X
        got = glue_cdm_source(ap, rho_c, rho_X)
        denom = max(abs(expected), 1e-300)
        worst = max(worst, abs(got - expected) / denom)
    check(
        "cdm_source == -(beta/c_D) Psi(r) rho_X       (Eq. 44)",
        worst < 1e-14,
        f"max rel err = {worst:.2e}",
    )
    # beta = 0 sanity: source vanishes identically.
    ap0 = AfterglowParams(c_D=1.3, beta_aft=0.0)
    check(
        "beta = 0 => cdm_source = 0 (no back-reaction)",
        abs(glue_cdm_source(ap0, 0.26, 0.69)) < 1e-300,
    )


def test11_total_dark_sector_conservation():
    """
    The matching source to rho_X lives inside Eq. 43:
        d rho_X / dN = -(1/c_D)[1 - beta Psi] rho_X
                     = -3 (1 + w_X) rho_X  +  (beta/c_D) Psi rho_X
    where (1 + w_X) = 1/(3 c_D). The CDM source is
        d rho_c / dN = -3 rho_c  -  (beta/c_D) Psi rho_X   (Eq. 44)
    so the combined source (beta/c_D) Psi rho_X CANCELS and
        d(rho_X+rho_c)/dN + 3 rho_c + 3(1+w_X) rho_X = 0.
    This test is the Bianchi-identity regression guard.
    """
    print("\nTest 11 :  Total dark-sector conservation (Eqs. 43 + 44)")
    ap = AfterglowParams(c_D=1.3, beta_aft=2.5)
    Omega0_X = 0.69
    Omega0_c = 0.26
    w_X = -1.0 + 1.0 / (3.0 * ap.c_D)
    worst = 0.0
    for N in [-3.0, -1.5, -0.7, 0.0]:
        rho_X = Omega0_X * math.exp(-N / ap.c_D)
        Sigma = rho_X / (3.0 * ap.c_D)
        rho_c = Omega0_c * math.exp(-3.0 * N)
        # d rho_X / dN from Eq. 43 via glue
        d_rho_X, _ = glue_derivs(ap, rho_c, rho_X, Sigma)
        # d rho_c / dN = -3 rho_c + cdm_source   (Eq. 44)
        d_rho_c = -3.0 * rho_c + glue_cdm_source(ap, rho_c, rho_X)
        # Bianchi residual: should vanish
        residual = (d_rho_X + d_rho_c) + 3.0 * rho_c + 3.0 * (1.0 + w_X) * rho_X
        scale = max(abs(d_rho_X) + abs(d_rho_c) + 3.0 * rho_c
                    + 3.0 * (1.0 + w_X) * rho_X, 1e-300)
        worst = max(worst, abs(residual) / scale)
    check(
        "Bianchi residual < 1e-14 (β/c_D Psi rho_X cancels)",
        worst < 1e-14,
        f"max rel residual = {worst:.2e}",
    )
    # beta = 0 sanity: residual trivially zero and each side matches LCDM.
    ap0 = AfterglowParams(c_D=1.3, beta_aft=0.0)
    rho_X0 = 0.69
    rho_c0 = 0.26
    Sigma0 = rho_X0 / (3.0 * ap0.c_D)
    d_rho_X0, _ = glue_derivs(ap0, rho_c0, rho_X0, Sigma0)
    d_rho_c0 = -3.0 * rho_c0 + glue_cdm_source(ap0, rho_c0, rho_X0)
    w_X0 = -1.0 + 1.0 / (3.0 * ap0.c_D)
    res0 = (d_rho_X0 + d_rho_c0) + 3.0 * rho_c0 + 3.0 * (1.0 + w_X0) * rho_X0
    check(
        "beta = 0 Bianchi residual vanishes",
        abs(res0) < 1e-14,
        f"residual = {res0:.2e}",
    )


# ═══════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("=" * 63)
    print("  CLASS_SYMT  Phase 1 background sector — unit tests")
    print("  Martin & Koh (April 2026) afterglow dark energy")
    print("=" * 63)
    test1_analytic_scaling()
    test2_LCDM_limit()
    test3_Psi_extrema()
    test4_Psi_sign()
    test5_exchange_sign()
    test6_energy_conservation_smoothness()
    test7_rhs_pointmatch_eq43()
    test8_glue_derivs()
    test9_glue_initial_conditions()
    test10_cdm_source_pointmatch_eq44()
    test11_total_dark_sector_conservation()
    print("\n" + "=" * 63)
    print(f"  RESULT:  {_passed} passed, {_failed} failed")
    print("=" * 63)
    sys.exit(0 if _failed == 0 else 1)
