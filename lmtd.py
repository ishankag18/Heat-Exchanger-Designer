"""
lmtd.py
-------
LMTD (Log Mean Temperature Difference) method for shell-and-tube
heat exchanger design, including the Bowman-Mueller-Nagle F correction
factor for 1-2 shell-and-tube configurations, generalized to N shell
passes using Fakheri's series relation.

References (standard textbook correlations, no proprietary data):
    - Kern, D.Q., "Process Heat Transfer"
    - Incropera, F.P., "Fundamentals of Heat and Mass Transfer"
"""

import math


def lmtd_counterflow(T_hot_in, T_hot_out, T_cold_in, T_cold_out):
    """
    Compute the counter-current LMTD.

    Parameters
    ----------
    T_hot_in, T_hot_out : float
        Hot fluid inlet / outlet temperature (deg C or K, consistent units).
    T_cold_in, T_cold_out : float
        Cold fluid inlet / outlet temperature.

    Returns
    -------
    float : LMTD (same units as inputs)
    """
    dT1 = T_hot_in - T_cold_out
    dT2 = T_hot_out - T_cold_in

    if dT1 <= 0 or dT2 <= 0:
        raise ValueError(
            "Non-physical temperature approach: check that hot > cold "
            "at both ends of the exchanger."
        )

    if abs(dT1 - dT2) < 1e-9:
        return dT1  # avoid divide-by-zero when dT1 == dT2

    return (dT1 - dT2) / math.log(dT1 / dT2)


def _r_p_parameters(T_hot_in, T_hot_out, T_cold_in, T_cold_out):
    """Compute R (capacity ratio) and P (thermal effectiveness) parameters."""
    R = (T_hot_in - T_hot_out) / (T_cold_out - T_cold_in)
    P = (T_cold_out - T_cold_in) / (T_hot_in - T_cold_in)
    return R, P


def _f_factor_single_shell(R, P):
    """
    Bowman-Mueller-Nagle F correction factor for ONE shell pass
    (with 2, 4, 6... tube passes).
    """
    if P <= 0:
        return 1.0
    if P >= 1:
        raise ValueError("P must be < 1 for a valid F-factor calculation.")

    if abs(R - 1.0) < 1e-9:
        # Special case R = 1
        numerator = P * math.sqrt(2)
        denom_log_arg = (2 - P * (2 - math.sqrt(2))) / (2 - P * (2 + math.sqrt(2)))
        denominator = (1 - P) * math.log(denom_log_arg)
        return numerator / denominator

    sqrt_term = math.sqrt(R ** 2 + 1)
    num = sqrt_term * math.log((1 - P) / (1 - P * R))
    log_arg = (2 - P * (R + 1 - sqrt_term)) / (2 - P * (R + 1 + sqrt_term))
    den = (R - 1) * math.log(log_arg)
    return num / den


def f_correction_factor(T_hot_in, T_hot_out, T_cold_in, T_cold_out, n_shell_passes=1):
    """
    Generalized LMTD correction factor F for N shell passes in series
    (Fakheri's approach: split the overall P into an equivalent
    per-shell P1, then apply the single-shell-pass formula).

    Parameters
    ----------
    n_shell_passes : int
        Number of shell passes in series (N). N=1 is the classic 1-2 TEMA E shell.

    Returns
    -------
    float : F, dimensionless (0 < F <= 1)
    """
    R, P = _r_p_parameters(T_hot_in, T_hot_out, T_cold_in, T_cold_out)
    N = n_shell_passes

    if N == 1:
        P1 = P
    else:
        if abs(R - 1.0) < 1e-9:
            P1 = P / (N - P * (N - 1))
        else:
            ratio = ((1 - P * R) / (1 - P)) ** (1.0 / N)
            P1 = (1 - ratio) / (R - ratio)

    F = _f_factor_single_shell(R, P1)
    return max(0.0, min(F, 1.0))


def required_area(Q, U, F, lmtd):
    """
    Solve for heat-transfer area from the design equation Q = U * A * F * LMTD.

    Parameters
    ----------
    Q : float   -> Duty (W)
    U : float   -> Overall heat transfer coefficient (W/m^2.K)
    F : float   -> LMTD correction factor
    lmtd : float -> Log mean temperature difference (K)

    Returns
    -------
    float : Required heat transfer area (m^2)
    """
    return Q / (U * F * lmtd)


def number_of_tubes(area, tube_od, tube_length):
    """
    Estimate number of tubes required for a given area.

    Parameters
    ----------
    area : float       -> total required area (m^2)
    tube_od : float    -> tube outer diameter (m)
    tube_length : float -> effective tube length (m)

    Returns
    -------
    int : number of tubes (rounded up)
    """
    area_per_tube = math.pi * tube_od * tube_length
    return math.ceil(area / area_per_tube)
