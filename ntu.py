"""
ntu.py
------
NTU-Effectiveness method for shell-and-tube heat exchangers.
Avoids iterative hand calculation by directly solving explicit
effectiveness correlations (Incropera et al.).
"""

import math


def heat_capacity_rates(m_hot, cp_hot, m_cold, cp_cold):
    """
    Compute hot/cold stream heat capacity rates and derived quantities.

    Parameters
    ----------
    m_hot, m_cold : float  -> mass flow rates (kg/s)
    cp_hot, cp_cold : float -> specific heats (J/kg.K)

    Returns
    -------
    dict with C_hot, C_cold, C_min, C_max, Cr
    """
    C_hot = m_hot * cp_hot
    C_cold = m_cold * cp_cold
    C_min = min(C_hot, C_cold)
    C_max = max(C_hot, C_cold)
    Cr = C_min / C_max
    return {"C_hot": C_hot, "C_cold": C_cold, "C_min": C_min, "C_max": C_max, "Cr": Cr}


def ntu_from_UA(U, A, C_min):
    """NTU = UA / C_min."""
    return (U * A) / C_min


def effectiveness_shell_and_tube(NTU, Cr, n_shell_passes=1):
    """
    Effectiveness for a shell-and-tube exchanger with N shell passes
    and 2N, 4N, ... tube passes (Incropera correlation).

    Parameters
    ----------
    NTU : float  -> overall number of transfer units (based on total UA)
    Cr : float   -> capacity ratio C_min/C_max (0 <= Cr <= 1)
    n_shell_passes : int -> number of shell passes N

    Returns
    -------
    float : effectiveness (0 to 1)
    """
    N = n_shell_passes
    NTU1 = NTU / N  # NTU per shell pass

    if NTU1 < 1e-9:
        return 0.0

    if Cr < 1e-9:
        # Condenser / evaporator limit
        eps1 = 1 - math.exp(-NTU1)
    else:
        sqrt_term = math.sqrt(1 + Cr ** 2)
        exp_term = math.exp(-NTU1 * sqrt_term)
        eps1 = 2.0 / (
            1 + Cr + sqrt_term * (1 + exp_term) / (1 - exp_term)
        )

    if N == 1:
        return eps1

    if Cr < 1e-9:
        return 1 - (1 - eps1) ** N

    if abs(Cr - 1.0) < 1e-9:
        return (N * eps1) / (1 + (N - 1) * eps1)

    ratio = ((1 - eps1 * Cr) / (1 - eps1)) ** N
    eps = (ratio - 1) / (ratio - Cr)
    return eps


def ntu_required_for_effectiveness(eps_target, Cr, n_shell_passes=1, tol=1e-8, max_iter=200):
    """
    Invert the effectiveness correlation to find the NTU required to hit
    a target effectiveness. Uses bisection (robust, no derivatives needed)
    since the closed-form correlation is only explicit in the NTU->eps
    direction.

    Returns
    -------
    float : NTU
    """
    lo, hi = 0.0, 50.0
    f_lo = effectiveness_shell_and_tube(lo, Cr, n_shell_passes) - eps_target
    f_hi = effectiveness_shell_and_tube(hi, Cr, n_shell_passes) - eps_target

    if f_hi < 0:
        raise ValueError(
            "Target effectiveness not achievable even at very high NTU "
            "for this Cr / shell-pass configuration."
        )

    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        f_mid = effectiveness_shell_and_tube(mid, Cr, n_shell_passes) - eps_target
        if abs(f_mid) < tol:
            return mid
        if f_mid > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def duty_from_effectiveness(eps, C_min, T_hot_in, T_cold_in):
    """Q = eps * C_min * (T_hot_in - T_cold_in)."""
    return eps * C_min * (T_hot_in - T_cold_in)


def outlet_temperatures(eps, C_min, C_max, T_hot_in, T_cold_in):
    """
    Recover outlet temperatures once effectiveness and Q are known.

    Returns
    -------
    dict with T_hot_out, T_cold_out, Q
    """
    Q = duty_from_effectiveness(eps, C_min, T_hot_in, T_cold_in)
    # Need actual C_hot, C_cold (not just min/max) passed in correct order.
    # This helper assumes caller supplies C_min/C_max already mapped to the
    # correct stream via the `hot_is_min` flag in HeatExchangerDesign.
    T_hot_out = T_hot_in - Q / C_min if C_min else T_hot_in
    T_cold_out = T_cold_in + Q / C_max if C_max else T_cold_in
    return {"T_hot_out": T_hot_out, "T_cold_out": T_cold_out, "Q": Q}
