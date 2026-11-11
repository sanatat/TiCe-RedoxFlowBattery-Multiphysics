# Research keyword coverage: University of Antwerp Ti–Ce project

This is a **traceability matrix**, not a claim of laboratory competency or quantitative experimental validation. The basis is Fu et al. (*ChemSusChem*, 2026, DOI [10.1002/cssc.202502732](https://doi.org/10.1002/cssc.202502732)) and its 41-page Supporting Information (SI). University of Antwerp / Hereijgers is the **application context**, not a collaborating or endorsing party.

**Three levels:** `Numerically implemented`, `Conditional scenario only`, and `Not established`.

| Target keyword / research gap | Files and reproducible evidence | Evidence level | Still missing for validated physical conclusions |
|---|---|---|---|
| **Mass Transfer** / **Convection–Diffusion–Reaction** | `transport_fv.py`, `porous_felt.py`, `extended_outputs/02_coupled_reactor_rate_flow_diffusion.csv`, `extensions/reactor.py` | Numerically implemented | True felt microgeometry, independent measured flux/limiting current |
| **Reactor Engineering** / **Darcy Flow** | `extensions/reactor.py`, `run_integrated.py`, `design_tradeoffs.py` | Numerically implemented as reduced order | Cell CAD, measured pressure drop, flow penetration, distributed current model |
| **Electrochemical Experiments** / **CV/EIS** | `cv_virtual.py`, `electrochem_cell.py`, SI S4/S5 table regression | Computational only | Physical potentiostat/galvanostat, CV, EIS acquisition, electrode assembly |
| **Ti–Ce Electrolyte Chemistry** | `electrolyte_chem.py`, `extensions/chemistry.py` | SI-referenced with conditional acid scenarios | Real metal–sulfate/MSA stability constants, ionic activity and redox speciation |
| **Mixed-Acid Electrolyte Formulation** | `extended_outputs/01_mixed_acid_conditional_screen.csv` | Conditional scenario only | Measured H2SO4/MSA effects and Ti/Ce redox/solubility behavior |
| **Electrolyte Precipitation / Solubility** | `extensions/chemistry.py` hypothetical conditional loading-capacity ratio | **Not established:** stress-test only | Chemical potentials, phase equilibria, precipitation kinetics, valid PHREEQC database |
| **Membrane Selection / Crossover** | `extensions/membrane.py`, `extended_outputs/03_membrane_fick_migration_water_scenarios.csv` | Conservative numerical model + scenarios | Ti permeability, membrane potential, complex charges, measured membrane A/L geometry, calibrated water flux |
| **Battery Cycling / Energy Efficiency** | `electrochem_cell.py`, original synthetic profiles | Numerically implemented | Raw voltage-time series, parameter discrimination, measured cells |
| **CFD / Flow Fields** | `lbm_flow.py`, `transport_fv.py`, scenario fields | Educational 2D numerical model | Measured channel/CAD, real felt flow topology, validation against flow maps |
| **Porous Electrode Engineering** | `porous_felt.py`, `extensions/reactor.py` | Reduced-order with assumed surface/porosity | BET/electrochemical surface area, felt permeability, true reaction distribution |
| **Scale-up / Cell Stacking** | `extensions/stack.py`, `extended_outputs/04_stack_scaling_and_manifold_scenarios.csv` | Conditional network model | Manifold engineering, shunt currents, stack endplates, compression, measured pump losses |
| **Long-Term Stability / Regeneration** | `extensions/lifetime.py`, `extended_outputs/05_lifetime_ideal_remix_vs_no_remix.csv` | Scenario only | Fit ageing and recovery to genuine time-series experiments, electrode kinetics and state transitions |

## Scientific equations represented

1. **Axial reactant balance:** `dC_react/dx = - I/(F Q L)` under imposed uniform geometric current; each tank's active species is handled separately.
2. **Nernst + Butler–Volmer:** local formal potential from redox ratios; inverted symmetric Butler–Volmer model for activation overpotential.
3. **Film mass transfer:** `j_lim = F k_m C_react` and concentration loss based on `log(1 - j/j_lim)` (the model returns infeasible rather than extrapolating beyond this bound).
4. **Reactor hydraulics:** Darcy pressure loss `Δp = μ L u / K` with unmeasured permeability `K`; two-loop pump power included.
5. **Membrane flux:** Fick permeability plus **linearized** effective-charge electromigration term and prescribed water transfer. Each metal's total moles are conserved by equal-and-opposite intertank exchange.
6. **Stack:** `V_stack = Σ V_cell`, with a common hydraulic pressure and parallel branch conductances under an assumed linear Darcy network.
7. **Degradation:** explicit reversible imbalance and irreversible felt/material loss with a conditional remix trigger. No inference about real degradation constants is claimed.

## Why the extensions do not completely close all gaps

The SI gives a geometric area of **9 cm²**, **20 mL** per tank, **20 mL/min** per half-cell, **5 mm** felt thickness and **N115** membrane. SI Tables S6–S8 supply experimental proxies for transport, kinetics and electrolyte properties at **MSA-only** conditions; they do **not** identify mixed-acid precipitation, the internal active surface area, carbon felt permeability, membrane field or lifetime rate constants. A model cannot deduce these uniquely from two energy-efficiency points or one polarization trace.

**Do not advertise:** experimentally performed CV/EIS, measured Ti–Ce battery building, verified mixed-acid precipitation, validated membrane electrodiffusion, a proven stack or laboratory-confirmed 3,721-cycle predictions.

**Defensible description:** “Built and verified an independent Python research prototype for reduced-order Ti–Ce electrochemical reactor modelling, 2D species transport, conditional membrane and mixed-acid scenarios, stack scaling and degradation sensitivity; published experimental parameters and unknown assumptions are explicitly separated.”

## Reproduction

```bash
python -m pip install -r requirements.txt
python run_extensions.py
python -m pytest -q tests/
```

Additional plots: `extended_figures/`. CSV source/output and scenario labels: `extended_outputs/`.


## Addendum (2026-11-11): experimental Ce solubility evidence and medium-specific equilibrium
`extensions/evidence_thermodynamics.py` and `run_evidence_thermodynamics.py` add patent-reported oxidation-state-resolved Ce solubility limits for *five exact acid mixtures* and literature Ce/SO4 conditional log beta specifically at 25 C and 2 M NaClO4. **Neither source establishes equilibrium precipitation Ksp for Ti/Ce in MSA/H2SO4.** No pH/SOC/temperature-wide solid-phase prediction is claimed; see [`THERMODYNAMICS_EVIDENCE.md`](THERMODYNAMICS_EVIDENCE.md).


## Reader-friendly overview

For the narrative connecting the technical methods to flow-battery research questions, see the main [`README.md`](README.md). For an exact keyword-to-file map rooted in the public Antwerp vacancy, see [`RESEARCH_KEYWORD_MAP.md`](RESEARCH_KEYWORD_MAP.md). These maps distinguish coded numerical work from unperformed laboratory experiments.
