# Attribution

CLASS_SYMT is an **independent** codebase derived from the public CLASS
(Cosmic Linear Anisotropy Solving System) source tree at:

    https://github.com/lesgourg/class_public

The base tree in this repository was seeded from a snapshot of
`class_public` (master branch) on 2026-04-09. The upstream git history
was intentionally stripped so this project maintains a clean first-parent
timeline tied to the afterglow dark-energy research program. This is a
source-level derivation, not a GitHub fork.

All files originally authored by the CLASS team remain under the CLASS
license. See `LICENSE` for the upstream terms. Modifications and new files
(anything under `source/afterglow/`, `include/afterglow/`, and this file)
are © 2026 Thomas Martin & Ingyu Koh, MIT licensed.

## Upstream acknowledgements

- Julien Lesgourgues and the CLASS developer team for the reference
  Boltzmann solver.
- CLASS papers I–IV (arXiv:1104.2932, 1104.2933, 1104.2934, 1104.2935).

## Research program

CLASS_SYMT supports the theory developed in:

- Martin, T. & Koh, I. *Dark Energy as the Thermodynamic Afterglow of a
  Hidden Gauge-Theory Transition.* April 2026.

Earlier iterations of the same research program (Glassy Dynamics κ model,
Yang-Mills condensate) were implemented on a GitHub fork of class_public
at `lawdroid/class_public`, branch `feature/kappa-evolution`. Those results
(notably the 7-parameter MCMC that ruled out the κ mechanism) informed the
design choices documented in `README_AFTERGLOW.md`.
