# Research keywords mapped to actual evidence

**Purpose.** This table maps the technical language of the University of Antwerp's 2026 Ti–Ce redox-flow-battery doctoral vacancy and related ELCAT research project to this **independent** codebase. It is an audit of *software functionality*, not a claim of experimental competence, collaboration with ELCAT, or experimentally verified predictions.

**Sources for research themes:** [University of Antwerp job posting](https://www.uantwerpen.be/en/jobs/vacancies/academic-staff/?descr=Doctoral-scholarship-holder-redox-flow-batteries&q=4562), [Jonas Hereijgers' research projects](https://www.uantwerpen.be/en/staff/jonas-hereijgers/research/), and [Fu et al., 2026](https://doi.org/10.1002/cssc.202502732).

| Keyword from ad / connected research topic | Code / repeatable evidence | Accurate coverage |
|---|---|---|
| **Redox flow batteries for renewable energy storage** | `electrochem_cell.py`, `run_integrated.py`, `extensions/stack.py` | Simulated Ti–Ce flow-battery research case; no physical device |
| **Mass transfer** | `transport_fv.py`, `porous_felt.py`, `figures/01_2d_mass_transfer.png` | Solved numerically for idealized geometry; conservation tested |
| **Reactor engineering / electrochemical reactors** | `extensions/reactor.py`, `extended_outputs/02_coupled_reactor_rate_flow_diffusion.csv` | Axial reduced-order electrochemical reactor and hydraulic scenario |
| **Titanium–cerium redox chemistry** | `electrolyte_chem.py`, `electrochem_cell.py`, `si_reference.json` | Published redox-state framework, formal Nernst and tank balances, not fully validated solution speciation |
| **Electrolyte formulation / type of acid / mixed-acid electrolytes** | `extensions/chemistry.py`, `extended_outputs/01_mixed_acid_conditional_screen.csv` | Conditional MSA/H2SO4 sensitivities; no demonstrated acid optimization |
| **Charge/discharge rate / current density / operating conditions** | `electrochem_cell.py`, `extensions/reactor.py`, `outputs/11_battery_rate_vs_paper.csv` | Numerical current/flow/cycling sensitivity; model assumptions remain |
| **Electrolyte stability / state of charge / temperature** | `electrolyte_chem.py`, `extensions/lifetime.py`, `thermo_evidence/` | Partial evidence-based parameter sweeps and scenario checks; not experimental stability mapping |
| **Flow field design / hydrodynamics / CFD** | `lbm_flow.py`, `transport_fv.py`, `figures/05_CFD_velocity_and_transport.png` | 2D illustrative LBM and species transport; no measured flow field/CAD |
| **Electrode geometry / porous felt / mass transport limitation** | `porous_felt.py`, `extensions/reactor.py`, `figures/03_porous_felt_reaction.png` | Porous reaction–diffusion; real internal area/tortuosity unknown |
| **Membrane / ion cross-over** | `extensions/membrane.py`, `extended_outputs/03_membrane_fick_migration_water_scenarios.csv` | Fick-type transport + conditional field/water terms; partial SI Ce anchor |
| **Power density / polarization losses / pumping cost** | `extensions/reactor.py`, `design_tradeoffs.py`, `extended_figures/03_coupled_reactor_polarization.png` | Reduced-order electrochemical/hydraulic scenarios; no full-cell validation |
| **Electrolyte precipitation / solubility** | `extensions/evidence_thermodynamics.py`, `thermo_evidence/outputs/` | Literature-constrained Ce solubility evidence only; **no real Ti/Ce precipitation prediction** |
| **Scale-up / stacking** | `extensions/stack.py`, `extended_outputs/04_stack_scaling_and_manifold_scenarios.csv` | Illustrative electrical/hydraulic network; no tested stack |
| **Battery lifetime / long-term stability / electrolyte regeneration** | `extensions/lifetime.py`, `extended_figures/07_long_term_regeneration_scenario.png` | Conditional fade/rebalancing models; not calibrated to 3,721 measured cycles |
| **CV / EIS / electrochemical experiments** | `cv_virtual.py`, `electrochem_cell.py`, `outputs/06_CV_SI_regression.csv` | Simulated voltammetry and impedance plus analysis of published CV; **no instrument experiments performed** |
| **Electrode surface functionalization / electrocatalysts** | No experimental or catalyst modeling module | **Not implemented**; relevant future research, not a completed claim |

## Interpretation

Three confidence levels are maintained throughout the repository:

- **Implemented numerical method:** the named solver runs and produces outputs, with tests of numerical properties where available. This does **not** make its assumed physical parameters real.
- **Conditional scenario / literature screening:** the code explores possible dependencies using disclosed parameters or applies published data only under their source conditions.
- **Not implemented:** physical electrochemistry experiments, surface functionalization and electrocatalyst development, a validated precipitation thermodynamic model, true 3D felt flow-field reconstruction, and experimentally validated stack lifetime.

For complete input provenance and model validity limits, read [`DATA_PROVENANCE.md`](DATA_PROVENANCE.md) and [`THERMODYNAMICS_EVIDENCE.md`](THERMODYNAMICS_EVIDENCE.md). The original experiment belongs to Fu *et al.*; all novel computational results and uncertainties in this repository are separate.