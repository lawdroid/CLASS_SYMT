"""
afterglow_theory.py
-------------------
Custom Cobaya theory wrapper for the CLASS_SYMT afterglow dark-energy
model (Martin & Koh, April 2026).

Inherits from the stock `classy` wrapper and adds two knobs that the
stock wrapper doesn't know about: c_D and beta_aft. These get mapped
to the CLASS input parameters that the afterglow glue layer reads.

Usage in a Cobaya yaml:
    theory:
      afterglow_theory.AfterglowTheory:
        python_path: ./
        extra_args:
          afterglow_on: "yes"

Stage-1 requests only background observables (H(z), D_L, D_A, theta_*);
Stage-2 additionally requests Cls and P(k).
"""
from __future__ import annotations
from cobaya.theories.classy import classy
from cobaya.log import LoggedError


class AfterglowTheory(classy):
    """
    Extends the stock Cobaya `classy` wrapper with afterglow parameters:
      * c_D      — dimensionless memory time (sets w_X = -1 + 1/(3 c_D))
      * beta_aft — loading-kernel coupling strength (Eq. 43 / Eq. 44)
      * Omega_X  — closure density (computed from 1 - sum_others in initialize_with_params)
    """

    # Declare that this theory can provide derived parameters on top of
    # what classy already offers.
    def must_provide(self, **requirements):
        # Inherit everything classy provides (H(z), rs, sigma_8, Cls, P(k))
        super().must_provide(**requirements)

    def calculate(self, state, want_derived=True, **params_values_dict):
        """
        Map MCMC sample -> CLASS input dictionary, then delegate to
        the stock classy compute.
        """
        # --- map sampled params to CLASS keys ---
        params = dict(params_values_dict)  # shallow copy
        c_D      = params.pop("c_D")
        beta_aft = params.pop("beta_aft")

        # Closure: Omega_X = 1 - Omega_b - Omega_cdm - Omega_r (radiation
        # fixed by T_cmb). CLASS handles this if we set Omega_Lambda = 0
        # and let afterglow carry the budget.
        params["Omega_Lambda"] = 0.0
        params["Omega_afterglow"] = -1.0       # sentinel -> compute from closure

        # --- afterglow inputs (read by src/background.c via input.c) ---
        params["afterglow_on"] = "yes"
        params["afterglow_c_D"] = c_D
        params["afterglow_beta"] = beta_aft
        params["afterglow_Sigma_today"] = -1.0  # auto: Omega0_X / (3 c_D)

        # Delegate to stock classy. This runs CLASS, fills state with
        # all standard observables, and we add afterglow-specific
        # derived quantities below.
        super().calculate(state, want_derived=want_derived, **params)

        if want_derived and state.get("derived") is not None:
            w_X = -1.0 + 1.0 / (3.0 * c_D)
            state["derived"]["w_X"] = w_X
            # Omega_X from CLASS's background table
            cosmo = self.classy
            bg = cosmo.get_background()
            # Locate afterglow rho_X column if present
            if "(.)rho_afterglow" in bg:
                rho_X_today = bg["(.)rho_afterglow"][-1]
                rho_crit_today = bg["(.)rho_crit"][-1]
                state["derived"]["Omega_X"] = rho_X_today / rho_crit_today
            else:
                # Fall back to closure
                O_b = params.get("omega_b", 0.0224) / (params.get("H0", 67.4) / 100.0) ** 2
                O_c = params.get("omega_cdm", 0.12) / (params.get("H0", 67.4) / 100.0) ** 2
                O_r = 9.2e-5          # photons + massless neutrinos
                state["derived"]["Omega_X"] = 1.0 - O_b - O_c - O_r

    def initialize(self):
        """
        Called once at start of chain. Check that the CLASS build
        actually has afterglow support compiled in.
        """
        super().initialize()
        try:
            # Probe: try setting an afterglow-only parameter. If the C
            # code doesn't know about it, classy raises.
            self.classy.set({"afterglow_on": "no"})
            self.classy.set({"afterglow_c_D": 1.0})
        except Exception as e:
            raise LoggedError(
                self.log,
                "Your CLASS binary does NOT have afterglow support. "
                "Rebuild CLASS_SYMT with patches 1-12 applied and "
                "afterglow_on/afterglow_c_D/afterglow_beta exposed in "
                "input.c. Error: %s" % e
            )
        self.log.info("AfterglowTheory initialized successfully.")
