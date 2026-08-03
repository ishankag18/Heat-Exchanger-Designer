"""
heat_exchanger.py
-----------------
Top-level design class that ties together the LMTD and NTU-Effectiveness
methods for a shell-and-tube heat exchanger, and automates calculation of:
    - Heat duty (Q)
    - Required area (A)
    - Number of tubes (N)
    - Overall heat transfer coefficient (U) [can be supplied or solved for]

Designed around a crude-oil preheat-train style problem: heating a cold
crude stream from a low temperature up towards a target using a hot
process stream, across one or more shell-and-tube exchangers.
"""

from dataclasses import dataclass, field
from typing import Optional

import lmtd
import ntu


@dataclass
class StreamData:
    name: str
    m_flow: float      # kg/s
    cp: float           # J/kg.K
    T_in: float         # deg C
    T_out: Optional[float] = None  # deg C (may be unknown / solved for)


@dataclass
class HeatExchangerDesign:
    hot: StreamData
    cold: StreamData
    U: float                     # overall heat transfer coefficient, W/m^2.K
    n_shell_passes: int = 1
    tube_od: float = 0.019        # m (3/4" tube, common in refinery service)
    tube_length: float = 6.0      # m

    def __post_init__(self):
        self._capacity = ntu.heat_capacity_rates(
            self.hot.m_flow, self.hot.cp, self.cold.m_flow, self.cold.cp
        )
        self.hot_is_min = self._capacity["C_hot"] <= self._capacity["C_cold"]

    # ------------------------------------------------------------------
    # LMTD METHOD
    # ------------------------------------------------------------------
    def solve_lmtd_method(self, Q):
        """
        Given a target duty Q (W), and known inlet/outlet temperatures on
        both streams, compute LMTD, F, required area, and number of tubes.

        Requires hot.T_out and cold.T_out to already be set (rating problem)
        or computed beforehand (e.g. from an energy balance).
        """
        if self.hot.T_out is None or self.cold.T_out is None:
            raise ValueError(
                "hot.T_out and cold.T_out must be known to use the LMTD "
                "method directly. Use solve_ntu_method() first if only "
                "one outlet temperature is known."
            )

        LMTD = lmtd.lmtd_counterflow(
            self.hot.T_in, self.hot.T_out, self.cold.T_in, self.cold.T_out
        )
        F = lmtd.f_correction_factor(
            self.hot.T_in, self.hot.T_out, self.cold.T_in, self.cold.T_out,
            n_shell_passes=self.n_shell_passes,
        )
        A = lmtd.required_area(Q, self.U, F, LMTD)
        N_tubes = lmtd.number_of_tubes(A, self.tube_od, self.tube_length)

        return {
            "LMTD": LMTD,
            "F": F,
            "corrected_LMTD": F * LMTD,
            "area_m2": A,
            "n_tubes": N_tubes,
            "Q_W": Q,
        }

    # ------------------------------------------------------------------
    # NTU-EFFECTIVENESS METHOD
    # ------------------------------------------------------------------
    def solve_ntu_method(self, area=None, target_T_cold_out=None):
        """
        Two modes:
          1) Rating mode: supply `area` -> compute NTU, effectiveness,
             duty, and both outlet temperatures directly (no iteration).
          2) Design mode: supply `target_T_cold_out` -> back-calculate the
             required effectiveness, then the required NTU and area.

        Returns a dict with Cr, NTU, effectiveness, Q, area, outlet temps.
        """
        Cmin = self._capacity["C_min"]
        Cmax = self._capacity["C_max"]
        Cr = self._capacity["Cr"]

        if area is not None:
            NTU_val = ntu.ntu_from_UA(self.U, area, Cmin)
            eps = ntu.effectiveness_shell_and_tube(NTU_val, Cr, self.n_shell_passes)
            Q = ntu.duty_from_effectiveness(eps, Cmin, self.hot.T_in, self.cold.T_in)

        elif target_T_cold_out is not None:
            Q_max = Cmin * (self.hot.T_in - self.cold.T_in)
            Q_target = self.cold.m_flow * self.cold.cp * (
                target_T_cold_out - self.cold.T_in
            )
            eps = Q_target / Q_max
            NTU_val = ntu.ntu_required_for_effectiveness(
                eps, Cr, self.n_shell_passes
            )
            area = (NTU_val * Cmin) / self.U
            Q = Q_target
        else:
            raise ValueError("Provide either `area` or `target_T_cold_out`.")

        # Map C_min/C_max back onto hot/cold sides to get real outlet temps
        if self.hot_is_min:
            T_hot_out = self.hot.T_in - Q / self._capacity["C_hot"]
            T_cold_out = self.cold.T_in + Q / self._capacity["C_cold"]
        else:
            T_hot_out = self.hot.T_in - Q / self._capacity["C_hot"]
            T_cold_out = self.cold.T_in + Q / self._capacity["C_cold"]

        n_tubes = lmtd.number_of_tubes(area, self.tube_od, self.tube_length)

        return {
            "Cr": Cr,
            "NTU": NTU_val,
            "effectiveness": eps,
            "Q_W": Q,
            "area_m2": area,
            "n_tubes": n_tubes,
            "T_hot_out": T_hot_out,
            "T_cold_out": T_cold_out,
        }

    # ------------------------------------------------------------------
    # CROSS-CHECK: solve both methods and compare
    # ------------------------------------------------------------------
    def design_summary(self, target_T_cold_out):
        """
        Full design workflow for a preheat-train style problem:
          1. NTU method solves for required area given a target cold outlet.
          2. LMTD method is used to cross-check the same duty/area using
             the resulting outlet temperatures.
        """
        ntu_result = self.solve_ntu_method(target_T_cold_out=target_T_cold_out)

        # Populate outlet temps so the LMTD cross-check has what it needs
        self.hot.T_out = ntu_result["T_hot_out"]
        self.cold.T_out = ntu_result["T_cold_out"]

        lmtd_result = self.solve_lmtd_method(Q=ntu_result["Q_W"])

        return {"ntu_method": ntu_result, "lmtd_method": lmtd_result}
