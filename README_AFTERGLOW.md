# CLASS_SYMT — CLASS with Hidden SU(2) Yang-Mills Afterglow Dark Energy

This is an **independent** Boltzmann code based on the
[CLASS](https://github.com/lesgourg/class_public) (Cosmic Linear Anisotropy
Solving System) reference implementation. It is **not a fork**: the upstream
git history has been stripped and development proceeds on a clean first-parent
timeline tied to the Martin & Koh (April 2026) paper
*"Dark Energy as the Thermodynamic Afterglow of a Hidden Gauge-Theory Transition."*

## What this code adds on top of CLASS

1. **Hidden SU(2)_D Yang-Mills dark sector**
   - Pre-confinement dark radiation bath `ρ_dr,D` (ΔN_eff knob)
   - Confinement epoch `z_conf` with smooth transfer source `Q_c`
2. **Causal afterglow fluid** (post-confinement)
   - Energy density `ρ_X`, pressure `p_X = −ρ_X + Σ`
   - Memory variable `Σ ≥ 0` (Müller–Israel–Stewart bulk stress)
   - Loading–unloading law: `u·∇Σ = −(Θ / 3c_D)[1 − β Ψ(r)] Σ`  (Eq. 37)
   - Derived exchange: `Q = −(β / c_D) H Ψ(r) ρ_X`           (Eq. 42)
   - Signed kernel: `Ψ(r) = 4r(r−1) / (1+r)^3`,  `r = ρ_c / ρ_X`
3. **Independent perturbation DOFs**
   - `δ_X`, `θ_X`, and crucially `δΣ` as a dynamical variable
     (the stock `fld` module collapses Σ algebraically and *loses the memory*)

## New parameters (explanatory.ini)

```
Omega_dr_dark     0.0        # pre-confinement dark radiation today-equivalent
z_conf            1.0e5      # confinement redshift
c_D               1.0        # Deborah number (relaxation / Hubble)
beta_aft          2.5        # loading coupling (0 = no crossing)
Sigma_today       auto       # Σ_0 = ρ_X_0 / (3 c_D) if auto
cs2_X             1.0        # rest-frame sound speed (causal = 1)
```

## Implementation phases

| Phase | Scope                                          | Est. time |
|-------|------------------------------------------------|-----------|
| 1     | `background.c`: ρ_X, Σ, ρ_dr,D ODEs            | ~1 week   |
| 2     | `perturbations.c`: δ_X, θ_X, δΣ (IS closure)   | ~2 weeks  |
| 3     | MCMC: Planck + DESI DR2 + Pantheon+ (+ SH0ES)  | ~1–2 weeks|

This initial commit is the **Phase 0 scaffold**: placeholder
`source/afterglow/afterglow.c` and `include/afterglow/afterglow.h` plus the
CLASS base tree. Phase 1 will follow on branch `feature/afterglow-background`.

## Lessons carried forward from the earlier Glassy Dynamics branch

The previous `feature/kappa-evolution` branch on `lawdroid/class_public` modified
Newton's constant via a frozen stiffness κ(z). A full 7-parameter MCMC
(Mar 2026) ruled that mechanism out:

- `κ_c = 0.998 [0.997, 0.999]` (68%)
- `H_0 = 66.8 [65.9, 67.7]` — data prefers LCDM, not κ=1.176
- `χ²/dof = 1.160` vs LCDM 1.166

**Takeaway:** modifying H(z) via κ shifts the acoustic peaks and the other
six parameters cannot absorb the damage. The afterglow model deliberately
leaves pre-recombination physics untouched and routes all late-time
modification through the causal fluid Σ, keeping r_s safe.

## Citing

If you use CLASS_SYMT, please cite both:

- Lesgourgues, J. *The Cosmic Linear Anisotropy Solving System (CLASS) I: Overview.*
  arXiv:1104.2932 (upstream CLASS code)
- Martin, T. & Koh, I. *Dark Energy as the Thermodynamic Afterglow of a
  Hidden Gauge-Theory Transition.* (2026)

## License

CLASS_SYMT inherits the CLASS license (MIT) for all code originally derived
from `class_public`. New afterglow files are © 2026 Martin & Koh, MIT licensed.
See `LICENSE` and `ATTRIBUTION.md`.
