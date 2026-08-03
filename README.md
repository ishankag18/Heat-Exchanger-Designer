# Heat Exchanger Design for Crude Preheat Train

A Python tool for designing shell-and-tube heat exchangers for a crude
oil preheat train (30°C → 250°C), using **both** the LMTD and
**NTU-Effectiveness** methods.

## Features

- **LMTD method** with automated **F correction factor** calculation for
  multi-shell-pass (TEMA E) configurations, generalized to N shells in
  series using the Fakheri series relation.
- **NTU-Effectiveness method** solved with closed-form correlations
  (Incropera) — no iterative hand calculation required — to compute
  effectiveness (ε) and capacity ratio (Cr) directly.
- Automated calculation of **required heat transfer area (A)**, **number
  of tubes (N)**, and design cross-checks between both methods.
- Design plots for optimisation:
  - Temperature profile along the exchanger
  - ε vs NTU curves for varying Cr
  - Area vs U sensitivity curve

## Project structure

```
heat_exchanger_design/
├── lmtd.py            # LMTD + F correction factor
├── ntu.py             # NTU-effectiveness correlations
├── heat_exchanger.py  # HeatExchangerDesign class (ties both methods together)
├── plotting.py        # Temperature profile, ε-NTU, area-vs-U plots
├── main.py            # Worked example: crude preheat train (30°C -> 250°C)
├── requirements.txt
└── README.md
```

## Quick start

```bash
pip install -r requirements.txt
python main.py
```

This runs a worked example (crude oil feed heated from 30°C, hot process
stream in at 280°C) and:

1. Solves the **NTU method** in design mode (target cold outlet temperature
   specified, required area back-calculated).
2. Cross-checks the same duty with the **LMTD method** using the resulting
   outlet temperatures, including the multi-shell-pass F factor.
3. Saves three plots to `outputs/`:
   - `temperature_profile.png`
   - `effectiveness_vs_ntu.png`
   - `area_vs_U.png`

## Using it for your own case

```python
from heat_exchanger import HeatExchangerDesign, StreamData

hot = StreamData(name="Hot stream", m_flow=45.0, cp=2600.0, T_in=280.0)
cold = StreamData(name="Crude feed", m_flow=60.0, cp=2100.0, T_in=30.0)

hx = HeatExchangerDesign(
    hot=hot, cold=cold,
    U=450.0,            # overall heat transfer coefficient, W/m^2.K
    n_shell_passes=2,
    tube_od=0.019,       # m
    tube_length=6.0,     # m
)

result = hx.design_summary(target_T_cold_out=220.0)
print(result["ntu_method"])
print(result["lmtd_method"])
```

Or, if you already know all four terminal temperatures (rating problem),
call `hx.solve_lmtd_method(Q=...)` and `hx.solve_ntu_method(area=...)`
directly.

## Methods implemented

- **LMTD**: counter-current LMTD, Bowman–Mueller–Nagle F-factor for a
  single 1-2 TEMA E shell, generalized to N shells in series via
  Fakheri's relation, `A = Q / (U·F·LMTD)`.
- **NTU-Effectiveness**: standard shell-and-tube correlation
  (Incropera, Eq. 11.29/11.30 style) for 1 shell pass (2, 4, 6... tube
  passes), extended to N shells in series.

## License

MIT — feel free to reuse for coursework, hackathon, or portfolio projects.
