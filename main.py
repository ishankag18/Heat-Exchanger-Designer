"""
main.py
-------
Worked example matching the project scope: designing a shell-and-tube
heat exchanger for a crude-oil preheat train, heating crude from 30°C
towards 250°C using a hot process stream, via both the LMTD and
NTU-Effectiveness methods, then generating the design plots.

Run:
    python main.py
"""

import os

from heat_exchanger import HeatExchangerDesign, StreamData
import plotting


def main():
    out_dir = "outputs"
    os.makedirs(out_dir, exist_ok=True)

    # --- Stream definitions -------------------------------------------------
    # Hot side: process/product stream giving up heat, e.g. atmospheric
    # residue returning from the column at high temperature.
    hot = StreamData(name="Hot process stream", m_flow=45.0, cp=2600.0, T_in=280.0)

    # Cold side: crude oil feed to be preheated
    cold = StreamData(name="Crude oil feed", m_flow=60.0, cp=2100.0, T_in=30.0)

    # --- Exchanger definition -------------------------------------------------
    hx = HeatExchangerDesign(
        hot=hot,
        cold=cold,
        U=450.0,             # W/m^2.K, typical fouled crude/hot-oil service
        n_shell_passes=2,     # e.g. two 1-2 TEMA E shells in series
        tube_od=0.019,        # 3/4" OD tubes
        tube_length=6.0,      # m
    )

    # This single exchanger represents one unit in the preheat train (the
    # full train spans 30C -> 250C across several exchangers in series).
    # Target chosen to be thermodynamically achievable for this duty
    # (effectiveness must stay below the Cr/NTU-limited maximum).
    target_cold_outlet = 180.0

    # --- Solve: NTU method first (design mode), then LMTD cross-check --------
    result = hx.design_summary(target_T_cold_out=target_cold_outlet)

    ntu_res = result["ntu_method"]
    lmtd_res = result["lmtd_method"]

    print("=" * 60)
    print("NTU-EFFECTIVENESS METHOD")
    print("=" * 60)
    print(f"Cr (capacity ratio)        : {ntu_res['Cr']:.4f}")
    print(f"NTU                        : {ntu_res['NTU']:.4f}")
    print(f"Effectiveness (epsilon)    : {ntu_res['effectiveness']:.4f}")
    print(f"Duty, Q (MW)               : {ntu_res['Q_W'] / 1e6:.3f}")
    print(f"Required area (m^2)        : {ntu_res['area_m2']:.2f}")
    print(f"Number of tubes            : {ntu_res['n_tubes']}")
    print(f"Hot outlet temperature (C) : {ntu_res['T_hot_out']:.2f}")
    print(f"Cold outlet temperature (C): {ntu_res['T_cold_out']:.2f}")

    print()
    print("=" * 60)
    print("LMTD METHOD (cross-check using solved outlet temperatures)")
    print("=" * 60)
    print(f"LMTD (uncorrected), K      : {lmtd_res['LMTD']:.3f}")
    print(f"F correction factor        : {lmtd_res['F']:.4f}")
    print(f"Corrected LMTD, K          : {lmtd_res['corrected_LMTD']:.3f}")
    print(f"Required area (m^2)        : {lmtd_res['area_m2']:.2f}")
    print(f"Number of tubes            : {lmtd_res['n_tubes']}")

    # --- Plots -----------------------------------------------------------
    plotting.plot_temperature_profile(
        T_hot_in=hot.T_in, T_hot_out=hot.T_out,
        T_cold_in=cold.T_in, T_cold_out=cold.T_out,
        save_path=os.path.join(out_dir, "temperature_profile.png"),
    )

    plotting.plot_effectiveness_ntu(
        n_shell_passes=hx.n_shell_passes,
        save_path=os.path.join(out_dir, "effectiveness_vs_ntu.png"),
    )

    plotting.plot_area_vs_U(
        Q=ntu_res["Q_W"], F=lmtd_res["F"], lmtd_value=lmtd_res["LMTD"],
        save_path=os.path.join(out_dir, "area_vs_U.png"),
    )

    print(f"\nPlots saved to ./{out_dir}/")


if __name__ == "__main__":
    main()
