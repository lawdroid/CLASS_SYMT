CLASS_SYMT: CLASS + Hidden SU(2)_D Yang-Mills Afterglow Dark Energy
===================================================================

CLASS_SYMT is an **independent** derivative of CLASS (not a fork) that
adds the hidden-sector "afterglow" dark-energy module of
Martin & Koh (April 2026),
**"Dark Energy as the Thermodynamic Afterglow of a Hidden Gauge-Theory
Transition."**

📄 **Reference paper (included in this repo):**
[`docs/Martin_Koh_Afterglow_April2026.pdf`](docs/Martin_Koh_Afterglow_April2026.pdf)

GitHub renders the PDF inline when you click the link above. The paper
is the single source of truth for every equation implemented in
`source/afterglow/` and `include/afterglow/`. Each comment in that code
cites the paper equation number it implements (all verified against the
PDF). See also:

- [`README_AFTERGLOW.md`](README_AFTERGLOW.md) — project overview, parameters, phase plan
- [`PLAN_AFTERGLOW.md`](PLAN_AFTERGLOW.md) — equation-to-code map and phased rollout
- [`ATTRIBUTION.md`](ATTRIBUTION.md) — upstream CLASS acknowledgements and license
- [`afterglow.ini`](afterglow.ini) — sample input file (LCDM-compatible by default)
- [`test/test_afterglow_bg.py`](test/test_afterglow_bg.py) — Phase 1a+1b unit tests (24 passing)
- [`PHASE1B_PATCHES.md`](PHASE1B_PATCHES.md) — 6 surgical patches for `background.c` / `input.c`
- [`.claude/SESSION_STATE.md`](.claude/SESSION_STATE.md) — crash-recovery handoff for continued sessions

**Interactive learning artifacts (open in any browser):**
- [`afterglow_simple.html`](afterglow_simple.html) — start here. Four quantities explained plainly: ρ_X, p_X, w, a. Slide c_D and N, watch ρ_X(a) = 0.70·a^(−1/c_D) evolve.
- [`afterglow_interactive_learn.html`](afterglow_interactive_learn.html) — deeper: Ψ(r) kernel, w_X(c_D), ρ evolution, Q exchange — Play / Learn / Quiz tabs.
- [`afterglow_class_connection.html`](afterglow_class_connection.html) — why `dρ_X/dloga` is THE connection to CLASS. Click-a-quantity + 5-step walk-through of one ndf15 integration step.

**Physics notation notes (carried over from 2026-04-14 Q&A):**
- Overdot `˙` = `d/dt` by convention. The glue function `afterglow_glue_derivs()` returns `d/d(log a) = d/dN` (no H factor on the RHS). Related by `d/dt = H · d/dN`. CLASS steps in log(a), so the d/dN form is what the integrator consumes.
- `c_D` — canonical benchmark **1**, data-preferred **[0.6, 2.0]**, hard floor **> 1/3** (else w_X < −1 = phantom, forbidden). Three views: confinement parameter / EOS dial (`w_X = −1 + 1/(3c_D)`) / dilution rate (`ρ_X ∝ a^(−1/c_D)`).
- `Ω_X` sets the **amount** of dark energy (≈ 0.70 today). `c_D` sets the **character** (stiffness → EOS). Two independent knobs.

Key paper equations implemented (all VERIFIED):

| Eq. | Section | Meaning |
|-----|---------|---------|
| (31) | §4.1 | Intrinsic EoS `w_X = −1 + 1/(3 c_D)` |
| (32) | §4.2 | Density ratio `r = ρ_c / ρ_X` |
| (35) | §4.2 | Signed kernel `Ψ(r) = 4 r (r−1) / (1+r)³` |
| (36) | §4.2 | Sign convention of Ψ |
| (37) | §4.3 | Covariant loading–unloading law |
| (38) | §4.3 | `Σ̇ = −(H/c_D)[1 − β Ψ(r)] Σ` |
| (42) | §4.3 | Exchange law `Q = −(β/c_D) H Ψ(r) ρ_X` |
| (43) | §4.3 | `ρ̇_X = −(H/c_D)[1 − β Ψ(r)] ρ_X` |
| (55) | §5.1 | Late-time scaling `ρ_X ∝ a^(−1/c_D)` |

---

Below is the original CLASS README, kept verbatim for attribution.

---

CLASS: Cosmic Linear Anisotropy Solving System  {#mainpage}
==============================================

Authors: Julien Lesgourgues, Thomas Tram, Nils Schoeneberg

with several major inputs from other people, especially Benjamin
Audren, Simon Prunet, Jesus Torrado, Miguel Zumalacarregui, Francesco
Montanari, Deanna Hooper, Samuel Brieden, Daniel Meinert, Matteo Lucca, etc.

For download and information, see http://class-code.net


Compiling CLASS and getting started
-----------------------------------

(the information below can also be found on the webpage, just below
the download button)

Download the code from the webpage and unpack the archive (tar -zxvf
class_vx.y.z.tar.gz), or clone it from
https://github.com/lesgourg/class_public. Go to the class directory
(cd class/ or class_public/ or class_vx.y.z/) and compile (make clean;
make class). You can usually speed up compilation with the option -j:
make -j class. If the first compilation attempt fails, you may need to
open the Makefile and adapt the name of the compiler (default: gcc),
of the optimization flag (default: -O4 -ffast-math) and of the OpenMP
flag (default: -fopenmp; this flag is facultative, you are free to
compile without OpenMP if you don't want parallel execution; note that
you need the version 4.2 or higher of gcc to be able to compile with
-fopenmp). Many more details on the CLASS compilation are given on the
wiki page

https://github.com/lesgourg/class_public/wiki/Installation

(in particular, for compiling on Mac >= 10.9 despite of the clang
incompatibility with OpenMP).

To check that the code runs, type:

    ./class explanatory.ini

The explanatory.ini file is THE reference input file, containing and
explaining the use of all possible input parameters. We recommend to
read it, to keep it unchanged (for future reference), and to create
for your own purposes some shorter input files, containing only the
input lines which are useful for you. Input files must have a *.ini
extension. We provide an example of an input file containing a
selection of the most used parameters, default.ini, that you may use as a
starting point.

If you want to play with the precision/speed of the code, you can use
one of the provided precision files (e.g. cl_permille.pre) or modify
one of them, and run with two input files, for instance:

    ./class test.ini cl_permille.pre

The files *.pre are suppposed to specify the precision parameters for
which you don't want to keep default values. If you find it more
convenient, you can pass these precision parameter values in your *.ini
file instead of an additional *.pre file.

The automatically-generated documentation is located in

    doc/manual/html/index.html
    doc/manual/CLASS_manual.pdf

On top of that, if you wish to modify the code, you will find lots of
comments directly in the files.

Python
------

To use CLASS from python, or ipython notebooks, or from the Monte
Python parameter extraction code, you need to compile not only the
code, but also its python wrapper. This can be done by typing just
'make' instead of 'make class' (or for speeding up: 'make -j'). More
details on the wrapper and its compilation are found on the wiki page

https://github.com/lesgourg/class_public/wiki

Plotting utility
----------------

Since version 2.3, the package includes an improved plotting script
called CPU.py (Class Plotting Utility), written by Benjamin Audren and
Jesus Torrado. It can plot the Cl's, the P(k) or any other CLASS
output, for one or several models, as well as their ratio or percentage
difference. The syntax and list of available options is obtained by
typing 'pyhton CPU.py -h'. There is a similar script for MATLAB,
written by Thomas Tram. To use it, once in MATLAB, type 'help
plot_CLASS_output.m'

Developing the code
--------------------

If you want to develop the code, we suggest that you download it from
the github webpage

https://github.com/lesgourg/class_public

rather than from class-code.net. Then you will enjoy all the feature
of git repositories. You can even develop your own branch and get it
merged to the public distribution. For related instructions, check

https://github.com/lesgourg/class_public/wiki/Public-Contributing

Using the code
--------------

You can use CLASS freely, provided that in your publications, you cite
at least the paper `CLASS II: Approximation schemes <http://arxiv.org/abs/1104.2933>`. Feel free to cite more CLASS papers!

Support
-------

To get support, please open a new issue on the

https://github.com/lesgourg/class_public

webpage!
