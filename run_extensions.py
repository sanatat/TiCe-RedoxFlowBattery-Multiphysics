"""Reproduce five documented, conditional research-gap extensions.

All plots are GENERATED here, not copied from Fu et al. (2026).  The mixed
acid/solubility, complex charges, felt properties, stack manifolds and ageing
rates are hypothetical; the reported SI values anchor ONLY specified inputs.
"""
from dataclasses import replace
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from extensions.chemistry import MixedAcidHypothesis,mixed_acid_scenario
from extensions.reactor import ReactorScenario,simulate_discharge
from extensions.membrane import MembraneScenario,simulate_membrane
from extensions.stack import StackScenario,stack_operating_point
from extensions.lifetime import LifetimeScenario,simulate_lifetime

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'extended_outputs';FIG=ROOT/'extended_figures'


def write_csv(name,rows):
    if not rows:return
    with (OUT/name).open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def newfig(x,y,title,xlabel,ylabel,filename,yscale=None):
    fig,ax=plt.subplots(figsize=(8.2,4.8))
    if isinstance(y,dict):
        for k,v in y.items():ax.plot(x,v,linewidth=2,label=k)
        ax.legend(loc='best',fontsize=8)
    else:ax.plot(x,y,linewidth=2)
    ax.set(title=title,xlabel=xlabel,ylabel=ylabel)
    if yscale:ax.set_yscale(yscale)
    ax.grid(True,alpha=.25)
    fig.tight_layout();fig.savefig(FIG/filename,dpi=170);plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True);FIG.mkdir(exist_ok=True)
    metadata={'source':'Fu et al. ChemSusChem 2026 DOI 10.1002/cssc.202502732, SI',
              'status':'Research prototype; unmeasured scenario parameters, no lab validation',
              'modules':{}}

    # 1. Mixed-acid and CONDITIONAL solubility screen.  These are not measurements
    # and do not use validated sulfate/methanesulfonate formation constants.
    acid_rows=[];fvals=np.linspace(0,1,21)
    for hyp_name,sigma_scale,mu_scale,beta in [
        ('illustrative_A',1.15,.9,1.),('illustrative_B',.75,1.35,.15),
        ('illustrative_C',1.45,.8,5.)]:
        for f in fvals:
            h=MixedAcidHypothesis(sulfate_fraction=float(f),
                sulfate_to_msa_conductivity_factor=sigma_scale,
                sulfate_to_msa_viscosity_factor=mu_scale,
                conditional_beta_Ce_M_inv=beta,
                conditional_beta_Ti_M_inv=beta*2)
            r=mixed_acid_scenario(h)
            acid_rows.append({'scenario':hyp_name,**{k:v for k,v in r.items() if isinstance(v,(int,float,bool))}})
    write_csv('01_mixed_acid_conditional_screen.csv',acid_rows)
    newfig(fvals,{name:[r['conductivity_mS_cm_scenario'] for r in acid_rows if r['scenario']==name]
                  for name in ('illustrative_A','illustrative_B','illustrative_C')},
           'Mixed-acid conductivity SENSITIVITY (not measured)','Nominal H2SO4 fraction','Conductivity (mS/cm)',
           '01_mixed_acid_conditional_conductivity.png')
    newfig(fvals,{name:[r['Ce_loading_to_assumed_capacity_ratio'] for r in acid_rows if r['scenario']==name]
                  for name in ('illustrative_A','illustrative_B','illustrative_C')},
           'Conditional solubility envelope — NOT precipitation prediction',
           'Nominal H2SO4 fraction','Loading / hypothetical capacity (-)',
           '02_conditional_solubility_NOT_prediction.png')
    metadata['modules']['mixed_acid']='Measured 2 M MSA anchor; hypothetical sulfate perturbations and conditional binding.'

    # 2. Axially coupled galvanostatic reactor + transport / rate sensitivity.
    reac=[]
    for ds in (1.0,.30,.10):
        for q in (10.,20.,30.):
            for j in (10.,20.,30.,40.,50.,60.,70.,80.,90.):
                s=ReactorScenario(apparent_D_scale=ds,flow_mL_min=q,j_mA_cm2=j)
                r=simulate_discharge(s)
                reac.append({'D_scale_of_SI_apparent':ds,'flow_mL_min':q,'j_mA_cm2':j,
                             'feasible':r['feasible'],
                             'cell_voltage_V':r.get('cell_voltage_V',np.nan),
                             'activation_loss_V':r.get('mean_activation_V',np.nan),
                             'mass_transfer_loss_V':r.get('mean_concentration_V',np.nan),
                             'ohmic_loss_V':r.get('mean_ohmic_V',np.nan),
                             'pressure_drop_Pa':r['pressure_drop_Pa'] if r['feasible'] else np.nan,
                             'pump_power_W':r['pump_power_W'],
                             'limiting_j_mA_cm2':r.get('minimum_local_limiting_j_mA_cm2',r.get('limiting_j_mA_cm2',np.nan))})
    write_csv('02_coupled_reactor_rate_flow_diffusion.csv',reac)
    jvals=np.array([10,20,30,40,50,60,70,80,90],float)
    newfig(jvals,{f'apparent D ×{ds:g}':[next(r['cell_voltage_V'] for r in reac if r['flow_mL_min']==20 and r['D_scale_of_SI_apparent']==ds and r['j_mA_cm2']==j) for j in jvals]
                 for ds in (1.,.30,.10)},
           'Axially coupled reactor: transport sensitivity (scenario)',
           'Applied current density (mA/cm2)','Predicted discharge voltage (V)',
           '03_coupled_reactor_polarization.png')
    nominal=simulate_discharge(ReactorScenario())
    newfig(nominal['x_m']*1000,
           {'Ce(IV) reactant':nominal['Ce4_profile_mol_m3'],
            'Ti(III) reactant':nominal['Ti3_profile_mol_m3']},
           'Coupled galvanostatic consumption along assumed 30-mm channel',
           'Channel coordinate (mm)','Local active concentrations (mol/m3)',
           '04_reactor_axial_concentration.png')
    metadata['modules']['reactor']='1D galvanostatic convection/consumption, local Nernst/BV/film/ohmic and Darcy loss; not full porous CFD.'
    metadata['reactor_nominal_50mA_cm2']={k:float(nominal[k]) for k in
         ('cell_voltage_V','mean_activation_V','mean_concentration_V','mean_ohmic_V',
          'pressure_drop_Pa','Ce_balance_relative_error','Ti_balance_relative_error')}

    # 3. Membrane Ce/Ti migration and optional water-flow scenario.
    membrane_rows=[];mem_traces={}
    mem_cases={
       'fick_only':MembraneScenario(),
       'charged_complex_z_minus_one':MembraneScenario(z_Ce_effective=-1.,potential_drop_V=.01),
       'water_transfer_assumed':MembraneScenario(water_transfer_m3_s=1e-12),
    }
    for label,s in mem_cases.items():
        r=simulate_membrane(s,hours=120,initial_concentrations_mol_m3=((500,0),(0,0)))
        mem_traces[label]=r['trace']
        for row in r['trace']:
            membrane_rows.append({'scenario':label,**row})
        assert abs(r['final_Ce_mol']-r['initial_Ce_mol'])<1e-12
    write_csv('03_membrane_fick_migration_water_scenarios.csv',membrane_rows)
    xx=[r['time_h'] for r in mem_traces['fick_only']]
    newfig(xx,{k:[row['Ce_B_M'] for row in v] for k,v in mem_traces.items()},
           'Ce crossover: same membrane P, unmeasured field/water scenarios',
           'Elapsed time (h)','Ce in receiver (mol/L)',
           '05_membrane_crossover_migration.png')
    metadata['modules']['membrane']='Conservative Fick + linearized effective-charge migration + prescribed water-volume transfer.'

    # 4. Stack series voltage, parallel hydraulics, conditional maldistribution.
    stack_rows=[]
    for n in (1,2,4,8,16):
        for resistance_spread in (0.,.2):
            s=StackScenario(n_cells=n,hydraulic_resistance_spread=resistance_spread)
            r=stack_operating_point(s)
            stack_rows.append({'n_cells':n,'hydraulic_resistance_spread':resistance_spread,
                'feasible':r['feasible'],'stack_voltage_V':r.get('stack_voltage_V',np.nan),
                'gross_power_W':r.get('gross_power_W',np.nan),
                'net_power_W':r.get('net_power_W',np.nan),
                'pump_power_W':r.get('pump_power_W',np.nan),
                'min_cell_voltage_V':r.get('min_cell_voltage_V',np.nan),
                'max_cell_voltage_V':r.get('max_cell_voltage_V',np.nan),
                'flow_balance_relative_error':r.get('flow_balance_relative_error',np.nan)})
    write_csv('04_stack_scaling_and_manifold_scenarios.csv',stack_rows)
    xx=[1,2,4,8,16]
    newfig(xx,{f'resistance spread {spread:.0%}':[next(r['net_power_W'] for r in stack_rows if r['n_cells']==n and r['hydraulic_resistance_spread']==spread) for n in xx]
               for spread in (0.,.2)},
           'Illustrative electrically series-connected Ti-Ce stack scaling',
           'Number of cells','Net electrical power after pump (W)',
           '06_stack_net_power.png')
    metadata['modules']['stack']='Electrically series, Darcy hydraulically parallel; manifold resistance spread assumed; no shunt currents.'

    # 5. Long-term irreversible loss vs reversible remix scenarios.
    life_rows=[]
    life_results={}
    for label,remix in (('no_regeneration',False),('threshold_remix',True)):
        lr=simulate_lifetime(LifetimeScenario(cycles=1200,remix_enabled=remix))
        life_results[label]=lr
        life_rows.extend({'scenario':label,**row} for row in lr['trace'])
    write_csv('05_lifetime_ideal_remix_vs_no_remix.csv',life_rows)
    xx=[r['cycle'] for r in life_results['no_regeneration']['trace']]
    newfig(xx,{k:[r['capacity_fraction'] for r in v['trace']] for k,v in life_results.items()},
           'Hypothetical long-term fade and capacity recovery — not paper data',
           'Simulated cycle number','Normalized discharge capacity (-)',
           '07_long_term_regeneration_scenario.png')
    metadata['modules']['lifetime']='Phenomenological reversible imbalance, irreversible felt/inventory fade; not fitted to 3721 cycles.'
    metadata['lifetime_remix_events_1200']=len(life_results['threshold_remix']['remix_cycles'])

    (OUT/'extension_run_summary.json').write_text(json.dumps(metadata,indent=2,ensure_ascii=False),encoding='utf-8')
    print('Wrote 5 CSV datasets, 7 generated charts, and run provenance:',OUT,FIG)
    print('Nominal cell voltage (scenario):',round(nominal['cell_voltage_V'],4),'V')
    print('All experimental-vs-hypothetical boundaries: MODEL_COVERAGE.md')
    return metadata


if __name__=='__main__':main()
