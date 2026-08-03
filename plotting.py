"""
plotting.py
-----------
Generates the design/optimisation plots referenced in the project:
    1. Temperature profile along the exchanger (hot & cold streams)
    2. Effectiveness (epsilon) vs NTU curves for different Cr values
    3. Area vs U sensitivity curve (for design optimisation)
"""

import numpy as np
import matplotlib.pyplot as plt

import ntu as ntu_module


def plot_temperature_profile(T_hot_in, T_hot_out, T_cold_in, T_cold_out,
                              n_points=100, save_path=None):
    """
    Plots a simplified linear temperature profile of hot and cold streams
    across the exchanger length (counter-current arrangement assumed).
    """
    x = np.linspace(0, 1, n_points)
    T_hot = T_hot_in + (T_hot_out - T_hot_in) * x
    T_cold = T_cold_out + (T_cold_in - T_cold_out) * x

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(x, T_hot, label="Hot stream", color="#c0392b", linewidth=2.2)
    ax.plot(x, T_cold, label="Cold stream (crude)", color="#2980b9", linewidth=2.2)
    ax.fill_between(x, T_hot, T_cold, color="grey", alpha=0.08)

    ax.set_xlabel("Normalized exchanger length (cold-inlet end = 0)")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title("Temperature Profile — Counter-Current Shell & Tube HX")
    ax.legend()
    ax.grid(alpha=0.3)

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def plot_effectiveness_ntu(Cr_values=(0.0, 0.25, 0.5, 0.75, 1.0),
                            n_shell_passes=1, ntu_max=5, save_path=None):
    """
    Classic epsilon-NTU chart for a shell-and-tube exchanger (N shell
    passes), one curve per Cr value.
    """
    ntu_range = np.linspace(0, ntu_max, 200)

    fig, ax = plt.subplots(figsize=(7, 5))
    for Cr in Cr_values:
        eps_vals = [
            ntu_module.effectiveness_shell_and_tube(n, Cr, n_shell_passes)
            for n in ntu_range
        ]
        ax.plot(ntu_range, eps_vals, label=f"Cr = {Cr:.2f}", linewidth=2)

    ax.set_xlabel("NTU")
    ax.set_ylabel("Effectiveness (ε)")
    ax.set_title(f"ε–NTU Curves (Shell-and-Tube, N = {n_shell_passes} shell pass)")
    ax.set_ylim(0, 1.02)
    ax.legend(title="Capacity ratio")
    ax.grid(alpha=0.3)

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def plot_area_vs_U(Q, F, lmtd_value, U_range=None, save_path=None):
    """
    Sensitivity plot: required area vs overall heat transfer coefficient U,
    useful for design optimisation / fouling-margin studies.
    """
    if U_range is None:
        U_range = np.linspace(150, 1200, 200)  # typical crude-oil service range, W/m2.K

    area_vals = Q / (U_range * F * lmtd_value)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(U_range, area_vals, color="#27ae60", linewidth=2.2)
    ax.set_xlabel("Overall Heat Transfer Coefficient, U (W/m²·K)")
    ax.set_ylabel("Required Area, A (m²)")
    ax.set_title("Design Sensitivity: Required Area vs U")
    ax.grid(alpha=0.3)

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig
