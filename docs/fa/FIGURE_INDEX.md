# Chart and data index

| Figure | Output | Covers |
|---|---|---|
| 01 | `01_2d_mass_transfer.png` | Species concentration field from conservative FV mass transfer |
| 02 | `02_mass_transfer_sweep.png` | Mass transfer vs flow and diffusion uncertainty |
| 03 | `03_porous_felt_reaction.png` | Reaction penetration into real-thickness 5-mm felt, assumed pore parameters |
| 04 | `04_pressure_vs_flow.png` | Reactor hydraulic cost vs flow and channel height |
| 05 | `05_CFD_velocity_and_transport.png` | Coupled LBM CFD velocity and FVM scalar transport with/without baffles |
| 06 | `06_published_CV_peak_regression.png` | Analysis of **published measured** SI CV peak tables S4/S5 |
| 07 | `07_virtual_CV.png` | Virtual CV from numerical diffusion/BV; not new experiments |
| 08 | `08_virtual_EIS.png` | Virtual small-signal equivalent-circuit EIS; not measured |
| 09 | `09_electrolyte_concentration_effect.png` | Measured SI S8 conductivity/viscosity |
| 10 | `10_chemistry_nernst_temperature.png` | Ideal Nernst SOC & measured SI S7 temperature effects |
| 11 | `11_membrane_crossover_sweep.png` | Crossover with SI Ce permeability and assumed H-cell geometry |
| 12 | `12_battery_cycling_and_ideal_remix.png` | 0D coupled cycling, high-P stress, ideal rebalancing |
| 13 | `13_battery_rate_comparison.png` | Model sensitivity vs two published EE anchors, not validated |

| 14 | `14_reported_DFT_coordination_energies.png` | Analysis of SI Table S2 existing DFT binding energies, **not** a new quantum calculation or solution speciation |


Detailed raw tables are in `outputs/*.csv` and structured summary in `outputs/summary.json`. Run `python run_integrated.py` to regenerate all files.
