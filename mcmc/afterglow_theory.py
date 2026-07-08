"""
afterglow_theory.py
-------------------
Custom Cobaya theory wrapper for the CLASS_SYMT afterglow dark-energy
model (Martin & Koh, April 2026).
"""
from __future__ import annotations
from cobaya.theories.classy import classy
from cobaya.log import LoggedError


class AfterglowTheory(classy):

    _custom_derived = (
        "w_X", "Omega_X",
        "sigma8_at_z0p1", "sigma8_at_z0p3", "sigma8_at_z0p5",
        "fsigma8_at_z0p1", "fsigma8_at_z0p3", "fsigma8_at_z0p5", "late_branch_guard",
    )

    def must_provide(self, **requirements):
        super().must_provide(**requirements)
        # Strip our custom names from derived_extra if classy added them
        if hasattr(self, "derived_extra"):
            self.derived_extra = [
                p for p in self.derived_extra if p not in self._custom_derived
            ]

    def calculate(self, state, want_derived=True, **params_values_dict):
        params = dict(params_values_dict)
        c_D      = params.pop("c_D")
        beta_aft = params.pop("beta_aft", 0.0)

        params["afterglow_on"] = 1
        params["c_D"]          = c_D
        params["beta_aft"]     = beta_aft
        params["Sigma_today"]  = -1.0

        # Cache before super so _get_derived_all (called from inside) sees it
        self._last_c_D = c_D
        self._last_beta_aft = beta_aft
        super().calculate(state, want_derived=want_derived, **params)

    def _get_derived_all(self, derived_requested=True):
        # Stash and strip custom names from BOTH lists so parent never asks
        # CLASS to resolve them.
        full_op = list(self.output_params)
        full_de = list(self.derived_extra) if hasattr(self, "derived_extra") else []
        try:
            self.output_params = [p for p in full_op if p not in self._custom_derived]
            if hasattr(self, "derived_extra"):
                self.derived_extra = [
                    p for p in full_de if p not in self._custom_derived
                ]
            std_derived, std_extra = super()._get_derived_all(
                derived_requested=derived_requested)
        finally:
            self.output_params = full_op
            if hasattr(self, "derived_extra"):
                self.derived_extra = full_de

        if not derived_requested:
            return std_derived, std_extra

        # Now compute custom derived from CLASS state
        cosmo = self.classy
        try:
            c_D = getattr(self, "_last_c_D", 1.0)
            std_derived["w_X"] = -1.0 + 1.0 / (3.0 * c_D)
            beta_v = getattr(self, "_last_beta_aft", 0.0)
            std_derived["late_branch_guard"] = c_D - (1.0 + 4.0 * beta_v) / 3.0
            h = cosmo.h()
            R8 = 8.0 / h
            for z, label in ((0.1, "0p1"), (0.3, "0p3"), (0.5, "0p5")):
                s8z = cosmo.sigma(R8, z)
                fz = cosmo.scale_independent_growth_factor_f(z)
                std_derived[f"sigma8_at_z{label}"] = float(s8z)
                std_derived[f"fsigma8_at_z{label}"] = float(fz * s8z)
            bg = cosmo.get_background()
            key = next(
                (k for k in bg if "rho_afterglow" in k or "rho_X" in k),
                None,
            )
            if key is not None and "(.)rho_crit" in bg:
                std_derived["Omega_X"] = float(
                    bg[key][-1] / bg["(.)rho_crit"][-1])
            else:
                std_derived["Omega_X"] = float("nan")
        except Exception as e:
            self.log.warning(
                "AfterglowTheory: custom-derived failed: %s", e)
            for name in self._custom_derived:
                std_derived.setdefault(name, float("nan"))

        return std_derived, std_extra

    def initialize(self):
        super().initialize()
        try:
            self.classy.set({"afterglow_on": 1})
            self.classy.set({"c_D": 1.0})
        except Exception as e:
            raise LoggedError(
                self.log,
                "Your CLASS binary does NOT have afterglow support. "
                "Rebuild CLASS_SYMT with afterglow source files. Error: %s" % e
            )
        self.log.info("AfterglowTheory initialized successfully.")
