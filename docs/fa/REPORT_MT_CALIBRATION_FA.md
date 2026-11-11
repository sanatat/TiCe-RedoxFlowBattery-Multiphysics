# Ti–Ce Redox Flow Battery — Mass Transfer Calibration and Identifiability Audit

**Reference:** Fu et al., *ChemSusChem* 2026, DOI: [10.1002/cssc.202502732](https://doi.org/10.1002/cssc.202502732). The authors' 41-page Supporting Information (not redistributed) is the quantitative reference for electrode/membrane/flow parameters.

**Important limitation:** This is a reproducible numerical study and **not a physically validated digital twin** or proof of experimental training. The Figure 6c readings were *manually approximated from a low-resolution plot*. The 0.50 M mixed-electrolyte discharge branch is used, **not** primary numeric instrument measurements.

## 1. What is directly grounded in the article and SI?

- SI experimental pp.2–3: 9 cm² geometric electrode area, 5 mm heat-treated carbon felt, 20 mL reservoir per side, 20 mL/min per side, N115 membrane, 2.0 M MSA, 1.8 V maximum charge and 0.8 V discharge stop.
- SI Table S6 (p.34): apparent diffusion coefficients and electrode reaction rate constants for Ti and Ce, extracted by the authors from porous-felt CV. **The values are not independent bulk diffusivity measurements.**
- SI Table S7/S8 (pp.35–36): electrolyte conductivity and viscosity with temperature and metal concentration.
- SI Fig.S20 (p.21): N115 area resistance before cycling = **0.0942 Ω cm²**, after cycling = 0.0785 Ω cm².
- Paper Fig.6c (p.9): 0.5 M mixed-electrolyte discharge polarization curve; used here as 14 visually estimated training/assessment points. **Assumed visual reading uncertainty ±0.04 V is an analysis choice, not reported experimental SD.**
- Paper text p.9: energy efficiency at 30 mA/cm² = **86.02%**, at 90 mA/cm² = **64.37%**. Those are separate full-cycle rate tests, **not** data used to fit the discharge-only curve.

## 2. What actually ran?

1. The **original project** still includes conservative 2D finite-volume advection–diffusion–wall reaction, 1D porous-felt reaction–diffusion (Thiele/effectiveness), simple 2D CFD, membrane mass balance, virtual CV/EIS, Nernst, and battery cycling.
2. The **new inverse layer** fits a reduced three-parameter discharge polarization model to 10 visually estimated 0.50 M mRFB points, reserving four for internal checks.
3. The new **FV transport sweep** solves Ti(III) and Ce(IV) at Q=10/20/40 mL/min and apparent D×{1,0.1,0.01}; all 18 solves are checked against global mass balance.
4. The **effective transport envelope** converts the fitted *curvature threshold* into a *conditional* effective k under assumed SOC values and compares it to the **illustrative** channel-only FV k. The ratio is **NOT** a measurement of active carbon felt surface area.
5. An **independent porous-felt scenario grid** varies specific area, unknown intrinsic k, and D to demonstrate why the same performance can arise from different microscopic parameters.
6. Conditional engineering scenarios scale the fitted curvature parameter with predicted k(Q), and report voltage/pump tradeoffs. **This extrapolation is an assumption**, because the paper does not supply flow-rate-dependent measured discharge curves.

## 3. Numerical results from this run

| Quantity | Result | Scientific interpretation |
|---|---:|---|
| Training points | 10 | Approximate manual readings, Figure 6c |
| Held out points | 4 | From the *same graph*, not independent experiment |
| Train RMSE | **0.0259 V** | Versus hand-digitized curve |
| Internal held-out RMSE | **0.0261 V** | Same curve and assumptions |
| Extrapolated intercept E* | **1.704 V** | Fitted lumped value, not direct OCV measurement |
| Fitted lumped polarization ASR | **3.934 Ω cm²** | Absorbs ohmic + kinetic + other effects; **not actual ohmic ASR** |
| Effective current-curvature scale j* | **275.9 mA/cm²** | **Not validated limiting current**; dominated by last visual datapoint |
| TI(III) channel-only k, D=reported ×1, Q=20 | **6.78×10⁻⁵ m/s** | Assumed flat 1-mm channel geometry, not real felt |
| Max mass balance residual, 18 FV cases | **~2.0×10⁻¹⁴** | Numerical conservation check, not measurement |
| N115 voltage drop at 90 mA/cm² | **0.00848 V** | From directly reported N115 ASR; excludes all other cell losses |
| Single-pass theoretical utilization at 90 mA/cm² | **5.0%** | 9cm², 20mL/min, 0.5M limiting active species; assumes ideal feed |

### The decisive non-identifiability check

When the visually estimated point at 275 mA/cm² is **excluded**, the fitted curvature parameter goes to the optimizer's **upper bound (100,000 mA/cm²)** and its predicted voltage at 275 mA/cm² rises to **~0.48 V** instead of the manually read **~0.30 V**. The calibration therefore **does not identify a physically meaningful limiting current** independent of uncertain high-current graph readout. The bootstrap interval for j* near the lower allowed boundary is **not** a real confidence interval for physical diffusion limitation.

### Why a Carbon Felt surface area cannot be uniquely calculated here

For a simplified porous slab, effective diffusivity `D_eff = D ε/τ` and `φ=L sqrt(a_s k_surface/D_eff)` produce `k_felt = sqrt(D_eff a_s k_surface)*tanh(φ)`; many combinations of specific internal area `a_s`, intrinsic surface kinetic constant, tortuosity, wetting and local SOC give the same effective current response. Further, the S6 kinetics came from CV of porous felt already; treating that coefficient as a separately measured intrinsic pore-surface rate while also fitting active surface area would double-count geometry.

**Needed to establish a physically identifiable estimate:** raw polarization points with recorded SOC, limiting-current scans at multiple flow rates and concentrations, electrode porosity/wetted electroactive area measurements, true felt flow-field geometry with pressure-drop vs Q, and (ideally) half-cell EIS/CV with iR correction. Then a joint charge-transfer and transport inversion can be attempted.

## 4. Relevance to the Antwerp PhD skills

| Skill | This project now demonstrates | Still not demonstrated |
|---|---|---|
| Mass Transfer | 2D conservative species transport and quantified k(D,Q), dimensionless diagnostics and parameter identifiability | Measurement-based actual pore-scale k of this felt |
| Reactor Engineering | Single-pass utilization, time/flow dependence, model pressure/pump tradeoffs | Measured Δp of real cell or manufacturing process engineering |
| Electrochemistry | Published CV analysis and virtual CV/EIS; discharge polarization fit | Actual potentiostat EIS/CV and battery fabrication |
| Ti–Ce Chemistry | Quantitative published MSA solution transport/spectator balance | Mixed-acid speciation and independently validated precipitation kinetics |
| Membrane & Crossover | SI H-cell P and measured membrane ASR in previous integrated code | Accurate geometry-resolved migration + water transfer |
| Cycling & Regeneration | Dynamic mass balances and virtual cycles in previous integrated code | Measured 3,721-cycle replication and degradation rates |
| CFD/Flow Fields | Numerical coupling of velocity and species equations | CAD and verified real flow-field pressure/concentration maps |

## 5. Files

- `calibrate_transport.py`: approximate Figure 6c inverse analysis, transport and felt sensitivity, summary and plots.
- `design_tradeoffs.py`: conditional flow/pumping scenarios. Do not treat flow predictions as validated measurements.
- `validation_outputs/`: all new CSV tables and JSON.
- `validation_figures/`: 8 new diagnostic figures.
- The publisher figure crop is intentionally excluded from the public repository. Verify visually at the cited DOI, main paper Figure 6c.
- `tests/`: previous 25 tests + 12 new checks (37 total).
- All previous scripts, plots and notebook are retained.

## 6. CV-safe phrasing

*Developed and tested a reproducible Python modeling extension for a published Ti–Ce redox-flow-battery study, including finite-volume transport sensitivity to flow rate and apparent diffusion, a porous-felt reaction–diffusion parameter study, digitized discharge-polarization analysis, internal holdout checks, and an explicit identifiability assessment.*

Do **not** claim: *calibrated true active area, validated actual cell CFD, independently measured diffusion, experimental CV/EIS, proven full-stack lifetime prediction*.
