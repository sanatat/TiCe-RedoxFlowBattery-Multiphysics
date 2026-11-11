"""Generate reproducible outputs for integrated Ti-Ce RFB training study.
Run: python run_integrated.py [--fast]
No network or proprietary software required.
"""
from pathlib import Path
from dataclasses import replace
from argparse import ArgumentParser
import csv,json,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from transport_fv import Channel,channel_transport,lbm_velocity
from electrochem_cell import (Battery,SI,SPECIES,F,mass_transport_inputs,
                              published_apparent_D,published_apparent_k,
                              electrolyte_props,channel_for_battery,
                              init_tanks,membrane_step,metal_inventory,
                              virtual_eis,run_cycles,ocv)
from porous_felt import Felt,effectiveness,finite_difference_flux
from cv_virtual import simulate as virtual_cv
from electrolyte_chem import S2,formal_cell_voltage_ideal_proton

HERE=Path(__file__).resolve().parent
PLOTS=HERE/'figures';DATA=HERE/'outputs'
PLOTS.mkdir(exist_ok=True);DATA.mkdir(exist_ok=True)

def fig(name):
    plt.tight_layout();plt.savefig(PLOTS/name,dpi=155,bbox_inches='tight');plt.close()

def csvout(name,rows):
    if not rows:return
    with open(DATA/name,'w',newline='') as fh:
        dw=csv.DictWriter(fh,fieldnames=list(rows[0]));dw.writeheader();dw.writerows(rows)

def compact(obj):
    if isinstance(obj,(np.integer,np.floating)):return obj.item()
    if isinstance(obj,(np.ndarray,list,tuple)):return [compact(v) for v in obj]
    if isinstance(obj,dict):return {k:compact(v) for k,v in obj.items()}
    return obj

def main(fast=False):
    b=Battery();ch=channel_for_battery(b)
    nx,ny=(70,26) if fast else (110,42)
    res={'provenance':'Fu et al. ChemSusChem 2026 DOI 10.1002/cssc.202502732; accompanying SI 41 pages',
         'scope':'Computational educational study; no new laboratory observations; no electrochemical EIS measurements or precipitation validation.',
         'parameter_origin':{'cell':'SI pp 2-3','D_k':'SI Table S6','viscosity_sigma':'SI S7-S8','P_Ce':'SI Fig S4',
                             'membrane_ASR':'SI Fig S20','felt_porosity_tortuosity_active_area':'UNMEASURED scenarios'}}
    reference={s:published_apparent_D(s) for s in SPECIES}
    km,rfields=mass_transport_inputs(b,nx=nx,ny=ny)
    res['transport_SI_D_m2_s']=reference
    res['transport_channel_k_m_s']={k:float(v) for k,v in km.items()}
    # 1) Channel transport: primary scientific missing gap
    t=rfields['Ti3'];cc=t['c']
    plt.figure(figsize=(9,4.2));ax=plt.gca()
    image=ax.pcolormesh(t['x_m']*1000,t['y_m']*1000,cc,shading='nearest')
    plt.colorbar(image,label='Ti(III) concentration (mol/m3)')
    ax.set(xlabel='Axial position x (mm)',ylabel='Height y (mm)',title='2-D conservative advection-diffusion, Ti(III) reaction at wall')
    fig('01_2d_mass_transfer.png')
    rows=[]
    for species,rr in rfields.items():
        rows.append({'species':species,'reported_apparent_D_m2_s':reference[species],
                     'k_wall_transfer_m_s':rr['k_apparent_m_s'],
                     'equivalent_j_lim_mA_cm2_at_0p25M':rr['j_equiv_mA_cm2'],
                     'inlet_mol_s':rr['in_mol_s'],'outlet_mol_s':rr['out_mol_s'],
                     'consumption_mol_s':rr['reacted_mol_s'],
                     'mass_balance_relative_error':rr['mass_balance_relative_error']})
    csvout('01_species_transport_balance.csv',rows)
    res['mass_balance_max_abs_relative_error']=max(abs(x['mass_balance_relative_error']) for x in rows)
    # 2) Mass transfer flow AND diffusion sensitivity
    flows=[5,10,20,40,80]
    transport_rows=[]
    for flow in flows:
        chx=Channel(flow_m3_s=flow*1e-6/60)
        for scale in [1,.1,.01]:
            result=channel_transport(chx,reference['Ti3']*scale,250.,nx=nx,ny=ny)
            transport_rows.append({'flow_mL_min':flow,'D_scale':scale,'D_m2_s':reference['Ti3']*scale,
                                   'k_m_s':result['k_apparent_m_s'],
                                   'molar_conversion_fraction':result['reacted_mol_s']/result['in_mol_s'],
                                   'j_equivalent_mA_cm2':result['j_equiv_mA_cm2'],
                                   'mass_error':result['mass_balance_relative_error']})
    csvout('02_mass_transfer_sensitivity.csv',transport_rows)
    plt.figure(figsize=(7.8,4.6))
    for scale in [1,.1,.01]:
        subset=[r for r in transport_rows if r['D_scale']==scale]
        plt.plot(flows,[x['k_m_s']*1e6 for x in subset],'-o',label=f'D = {scale:g} × published apparent value')
    plt.xlabel('Each-cell electrolyte flow (mL/min)');plt.ylabel('Apparent wall transfer coefficient (µm/s)')
    plt.title('Mass-transfer sensitivity to flow and uncertain diffusion');plt.legend(fontsize=8);fig('02_mass_transfer_sweep.png')
    # 3) Porous electrode reaction diffusion, analytic vs independent FVM
    felt_rows=[]
    plt.figure(figsize=(7.9,4.6))
    for area in [2500,25000,100000]:
        f=Felt(active_area_m2_m3=area)
        item=effectiveness('Ti3',f)
        verify=finite_difference_flux('Ti3',f,n=500)
        coupled=channel_transport(ch,reference['Ti3'],250.,k_wall_m_s=item['k_wall_m_s'],nx=nx,ny=ny)
        felt_rows.append({'active_area_m2_per_m3_UNMEASURED':area,'porosity_ASSUMED':f.porosity,
                          'tortuosity_ASSUMED':f.tortuosity,'felt_thickness_m_SI':f.thickness_m,
                          'Thiele':item['Thiele'],'effectiveness':item['effectiveness'],
                          'reaction_penetration_depth_um':item['penetration_depth_m']*1e6,
                          'k_effective_m_s':item['k_wall_m_s'],'FVM_analytic_relative_error':verify['relative_error'],
                          'coupled_channel_j_mA_cm2':coupled['j_equiv_mA_cm2'],
                          'coupled_balance_error':coupled['mass_balance_relative_error']})
        plt.plot(item['z_m']*1e3,item['concentration_mol_m3']/250,
                 label=f'Active area {area:,} m²/m³; effectiveness {item["effectiveness"]:.3f}')
    plt.xlabel('Depth into 5-mm carbon felt (mm)');plt.ylabel('Reactant concentration / inlet concentration')
    plt.title('1-D porous felt diffusion + reaction (assumed pore parameters)');plt.legend(fontsize=8);fig('03_porous_felt_reaction.png')
    csvout('03_porous_felt_reaction_sweep.csv',felt_rows)
    res['porous_felt_independent_FVM_max_error']=max(x['FVM_analytic_relative_error'] for x in felt_rows)
    # 4) Reactor hydraulics and residence
    reactors=[]
    for flow in flows:
        for height in [.0007,.001,.0015]:
            p=Channel(height_m=height,flow_m3_s=flow*1e-6/60)
            transport=channel_transport(p,reference['Ti3'],250.,nx=nx,ny=ny)
            reactors.append({'flow_mL_min':flow,'assumed_height_mm':height*1000,
                             'Re':p.reynolds,'pressure_drop_Pa':p.pressure_drop_Pa,
                             'pump_power_both_cells_W':p.pump_power_two_cells_W,
                             'channel_residence_s':p.length_m*p.height_m*p.width_m/p.flow_m3_s,
                             'reservoir_turnover_s':b.tank_m3/p.flow_m3_s,
                             'Ti3_conversion_fraction':transport['reacted_mol_s']/transport['in_mol_s'],
                             'Ti3_k_m_s':transport['k_apparent_m_s']})
    csvout('04_reactor_pressure_flow.csv',reactors)
    plt.figure(figsize=(8,4.5))
    for ht in [.7,1.,1.5]:
        sub=[r for r in reactors if abs(r['assumed_height_mm']-ht)<1e-7]
        plt.plot(flows,[r['pressure_drop_Pa'] for r in sub],'-o',label=f'{ht:g} mm assumed channel')
    plt.xlabel('Flow rate per cell (mL/min)');plt.ylabel('Idealized laminar pressure drop (Pa)')
    plt.title('Reactor engineering: more flow reduces residence but costs pump head');plt.legend();fig('04_pressure_vs_flow.png')
    res['reactor_nominal']={'Re':ch.reynolds,'pressure_drop_Pa_assumed':ch.pressure_drop_Pa,
                            'pump_power_two_cells_W_assumed':ch.pump_power_two_cells_W,
                            'reservoir_turnover_s':b.tank_m3/b.flow_m3_s}
    # 5) 2-D Lattice Boltzmann flow coupled to FVM scalar (scenario only)
    lobm=[];fields=[];lbsize=(75,27) if fast else (90,31)
    for baffle in [False,True]:
        u,v,solid,raw=lbm_velocity(ch,nx=lbsize[0],ny=lbsize[1],baffles=baffle,steps=1350 if fast else 1850)
        rr=channel_transport(ch,reference['Ti3'],250.,flow_field=(u,v,solid))
        fields.append((u,v,solid,rr))
        lobm.append({'baffle_scenario':baffle,'raw_inlet_flow_lattice':float(raw['ux'][:,0].sum()) if 'ux' in raw else None,
                     'scaled_flow_m3_s':float(np.sum(u[:,0])*(ch.height_m/u.shape[0])*ch.width_m),
                     'transport_k_m_s':rr['k_apparent_m_s'],'reacted_mol_s':rr['reacted_mol_s'],
                     'mass_balance_error':rr['mass_balance_relative_error']})
    resistance_ratio=lobm[0]['raw_inlet_flow_lattice']/lobm[1]['raw_inlet_flow_lattice']
    for row in lobm:row['baffle_hydraulic_resistance_relative_to_open_channel_ESTIMATE']=resistance_ratio if row['baffle_scenario'] else 1.
    csvout('05_CFD_LBM_geometry_scenarios.csv',lobm)
    fig_,axs=plt.subplots(2,2,figsize=(10,6.8))
    for i,(u,v,solid,rr) in enumerate(fields):
        a=axs[i,0];img=a.imshow(np.ma.masked_array(u,mask=solid),origin='lower',aspect='auto')
        fig_.colorbar(img,ax=a,label='u (m/s)');a.set_title('LBM velocity - '+('baffles' if i else 'empty'))
        a=axs[i,1];img=a.imshow(np.ma.masked_array(rr['c'],mask=solid),origin='lower',aspect='auto',vmin=0,vmax=250)
        fig_.colorbar(img,ax=a,label='C Ti(III) (mol/m³)');a.set_title('FVM concentration using LBM velocity')
    fig_.suptitle('Illustrative 2-D flow/transport coupling: periodic lattice -> mapped cell');fig('05_CFD_velocity_and_transport.png')
    res['CFD_scenarios']=lobm
    # 6) Published CV tables and separate synthetic CV
    cv_rows=[]
    fig_,axs=plt.subplots(1,2,figsize=(10,4.3))
    for k,sp in enumerate(['Ce','Ti']):
        src=SI['CV_'+sp]
        rates=np.asarray(src['scan_mV_s']);x=np.sqrt(rates/1000.);y=np.asarray(src['Ipa_mA'])
        fit=np.polyfit(x,y,1);r2=1-np.sum((np.polyval(fit,x)-y)**2)/np.sum((y-y.mean())**2)
        cv_rows.append({'ion':sp,'anodic_peak_slope_mA_sqrt_s_per_sqrtV':fit[0],
                        'anodic_peak_intercept_mA':fit[1],'R2_fit_on_SI_S4_or_S5':r2,
                        'peak_delta_E_at_10mV_s_V':src['Epa_V_SHE'][0]-src['Epc_V_SHE'][0]})
        axs[k].plot(x,y,'o',label=f'Published SI data (R²={r2:.3f})')
        xx=np.linspace(x.min(),x.max(),60);axs[k].plot(xx,np.polyval(fit,xx),label='Linear regression')
        axs[k].set(xlabel='√ scan rate (√V/s)',ylabel='Peak oxidation current (mA)',title=sp+' SI CV peak regression')
        axs[k].legend(fontsize=8)
    fig('06_published_CV_peak_regression.png');csvout('06_CV_SI_regression.csv',cv_rows)
    fig_,axs=plt.subplots(1,2,figsize=(10,4.2))
    for k,sp in enumerate(['Ce','Ti']):
        if sp=='Ce': D=published_apparent_D('Ce3');k0=published_apparent_k('Ce3')
        else:D=published_apparent_D('Ti3');k0=published_apparent_k('Ti3')
        dat=virtual_cv(D,k0,conc=100.,scan=.005,ny=180 if fast else 260)
        axs[k].plot(dat[:,0],dat[:,1]/10)
        axs[k].set(xlabel='Potential relative to formal E° (V)',ylabel='Current density (mA/cm²)',title=sp+' virtual CV (simplified diffusion + BV)')
        csvout(f'06_virtual_CV_{sp}.csv',[{'potential_relative_V':a,'current_mA_cm2':j/10} for a,j in dat])
    fig('07_virtual_CV.png')
    res['CV_from_SI']=cv_rows
    # 7) EIS synthetic, not from measured SI data
    eis=virtual_eis(b,km)
    z=eis['Z_area_ohm_m2']*1e4 # to ohm cm2
    fig_,axs=plt.subplots(1,2,figsize=(10,4.2))
    axs[0].plot(z.real,-z.imag);axs[0].set(xlabel="Z' (Ω cm²)",ylabel="−Z'' (Ω cm²)",title='Virtual EIS Nyquist (not measured)')
    axs[1].semilogx(eis['frequency_Hz'],abs(z));axs[1].set(xlabel='Frequency (Hz)',ylabel='|Z| (Ω cm²)',title='Virtual EIS Bode magnitude')
    fig('08_virtual_EIS.png')
    csvout('07_virtual_EIS.csv',[{'frequency_Hz':f,'Zreal_ohm_cm2':zi.real,'Minus_Zimag_ohm_cm2':-zi.imag}
                                  for f,zi in zip(eis['frequency_Hz'],z)])
    res['EIS_assumed_R0_ohm_cm2']=eis['Rohm_area_ohm_m2']*1e4
    # 8) electrolyte, measured temperature and concentration SI data
    cc=SI['conductivity_viscosity_vs_concentration'];temp=SI['conductivity_viscosity_vs_temperature']
    concentration_rows=[]
    for x,s,mu in zip(cc['M_each_metal'],cc['conductivity_mS_cm'],cc['viscosity_cP']):
        p=Channel(mu_Pa_s=mu*1e-3)
        concentration_rows.append({'Ce_and_Ti_M_each':x,'SI_conductivity_mS_cm':s,'SI_viscosity_cP':mu,
                                  'idealized_pump_power_W_at_Q20':p.pump_power_two_cells_W,
                                  'bulk_ohmic_ASR_ohm_cm2_ASSUMED_0p5mm':(.0005/(s*.1))*1e4})
    csvout('08_electrolyte_concentration_SI_S8.csv',concentration_rows)
    fig_,axs=plt.subplots(1,2,figsize=(10,4.2))
    axs[0].plot(cc['M_each_metal'],cc['conductivity_mS_cm'],'-o');axs[0].set(xlabel='Each metal (M)',ylabel='Conductivity (mS/cm)',title='Measured SI S8: ionic conductivity')
    axs[1].plot(cc['M_each_metal'],cc['viscosity_cP'],'-o');axs[1].set(xlabel='Each metal (M)',ylabel='Viscosity (cP)',title='Measured SI S8: viscosity')
    fig('09_electrolyte_concentration_effect.png')
    # Nernst voltage vs SOC, not a chemical speciation/precipitation model
    soc=np.linspace(.05,.95,91)
    fig_,axs=plt.subplots(1,2,figsize=(10,4.2))
    for T in [298.15,308.15,318.15]:
        bb=replace(b,temp_K=T)
        axs[0].plot(soc,[ocv(init_tanks(bb,s),bb) for s in soc],label=f'{T-273.15:.0f}°C')
    axs[0].set(xlabel='SOC (idealized same on each half)',ylabel='Open-circuit cell potential (V)',title='Formal Nernst model, not speciation')
    axs[0].legend()
    tt=[25,30,35,45,65]
    axs[1].plot(tt,temp['conductivity_mS_cm'],'-o',label='conductivity (mS/cm)')
    ax2=axs[1].twinx();ax2.plot(tt,temp['viscosity_cP'],'--s',label='viscosity (cP)')
    axs[1].set(xlabel='Temperature (°C); room assumed 25°C',ylabel='Conductivity (mS/cm)',title='Measured SI S7 properties')
    ax2.set_ylabel('Viscosity (cP)')
    fig('10_chemistry_nernst_temperature.png')
    # 8b) Coordination models from DFT S2 (reported computations, not rerun DFT)
    fig_,axs=plt.subplots(1,2,figsize=(10,4.4))
    coord=[]
    for name,dat in S2.items():
        target=axs[0] if name.startswith('Ce') else axs[1]
        target.plot(dat['ligands'],dat['binding_Ha'],'-o',label=name)
        for lig,energy in zip(dat['ligands'],dat['binding_Ha']):
            coord.append({'reported_species':name,'methanesulfonate_ligands_x':lig,
                          'reported_electronic_binding_energy_Ha':energy,
                          'source':'SI Table S2; not an equilibrium formation constant'})
    for ax in axs:
        ax.set(xlabel='Number of MSA sulfonate ligands x',ylabel='Reported electronic binding ΔE (Hartree)')
        ax.legend(fontsize=8)
    fig_.suptitle('Published SI Table S2: coordination DFT energies (NOT equilibrium speciation)')
    fig('14_reported_DFT_coordination_energies.png')
    csvout('08b_SI_DFT_coordination_binding.csv',coord)
    h=[.5,1.,2.,3.]
    csvout('08b_formal_acid_activity_sensitivity.csv',[
        {'relative_proton_activity_ASSUMED':x,
         'formal_voltage_V_at_SOC_0p5':formal_cell_voltage_ideal_proton(.5,298.15,x),
         'note':'ideal only; MSA analytical acid concentration is NOT proton activity'} for x in h])
    res['acid_activity_ideal_cell_V']={str(x):formal_cell_voltage_ideal_proton(.5,298.15,x) for x in h}
    # 9) Membrane crossover: Fick diffusion and mass balance
    concentrations=[]
    for factor in [.1,1,10]:
        bb=replace(b,membrane_ce_m_s=b.membrane_ce_m_s*factor,membrane_ti_m_s=2e-10)
        tanks=init_tanks(bb,soc=.8)
        before=metal_inventory(tanks,bb)
        for hour in range(0,121):
            if hour>0: membrane_step(tanks,bb,3600)
            concentrations.append({'P_Ce_relative':factor,'time_h':hour,
                                   'Ce4_on_negative_M':tanks[1,1]/1000,
                                   'Ce4_on_positive_M':tanks[0,1]/1000,
                                   'total_Ce_mol':metal_inventory(tanks,bb)['Ce_mol'],
                                   'total_Ti_mol':metal_inventory(tanks,bb)['Ti_mol']})
        after=metal_inventory(tanks,bb)
        assert all(abs(before[k]-after[k])<1e-11 for k in before)
    csvout('09_membrane_crossover.csv',concentrations)
    plt.figure(figsize=(8,4.5))
    for fact in [.1,1,10]:
        sr=[r for r in concentrations if r['P_Ce_relative']==fact]
        plt.plot([r['time_h'] for r in sr],[r['Ce4_on_negative_M'] for r in sr],label=f'{fact:g} × SI P_Ce')
    plt.xlabel('Hours');plt.ylabel('Ce(IV) receiver reservoir (M)')
    plt.title('Membrane crossover sensitivity: NOT experimental H-cell geometry')
    plt.legend();fig('11_membrane_crossover_sweep.png')
    # 10) Battery cycle coupling, stress tests and ideal regeneration
    curves={};cycle_rows=[]
    for name,settings,regen in [('reported_P_baseline',b,()),
                               ('high_crossover_stress',replace(b,membrane_ce_m_s=b.membrane_ce_m_s*300),()),
                               ('high_crossover_with_ideal_remix',replace(b,membrane_ce_m_s=b.membrane_ce_m_s*300),(4,8))]:
        track,summ,final,inv=run_cycles(settings,km,cycles=12,j_mA_cm2=50,regenerate_after=regen)
        csvout('10_cycle_trace_'+name+'.csv',track)
        for d in summ:cycle_rows.append({'scenario':name,**d})
        curves[name]=summ
        res['inventory_'+name]=inv
    csvout('10_battery_cycle_summaries.csv',cycle_rows)
    fig_,axs=plt.subplots(1,2,figsize=(10,4.2))
    for name,ss in curves.items():
        axs[0].plot([r['cycle'] for r in ss],[r['discharge_Ah']*1000 for r in ss],'-o',ms=3,label=name)
        axs[1].plot([r['cycle'] for r in ss],[r['EE_percent'] for r in ss],'-o',ms=3,label=name)
    axs[0].set(xlabel='Cycle',ylabel='Discharge capacity (mAh)',title='Scenario capacity retention')
    axs[1].set(xlabel='Cycle',ylabel='Electrical energy efficiency (%)',title='Scenario EE (model only)')
    axs[0].legend(fontsize=7);fig('12_battery_cycling_and_ideal_remix.png')
    # 11) Rate test: comparison with article anchors, but NOT calibration to those anchors
    rate=[]
    for multiplier in [1,2,3]:
        local_km={k:float(v)*multiplier for k,v in km.items()}
        for j in [30,50,70,90]:
            _,ss,_,_=run_cycles(b,local_km,cycles=2,j_mA_cm2=j)
            rate.append({'transport_area_multiplier_ASSUMED':multiplier,'j_mA_cm2':j,
                         'second_cycle_model_EE_percent':ss[-1]['EE_percent'],
                         'second_cycle_model_CE_percent':ss[-1]['CE_percent'],
                         'published_EE_percent_exact_rate_protocol':({30:86.02,90:64.37}).get(j),
                         'note':'Not same initial SOC / preconditioning / reaction area; not validated'})
    csvout('11_battery_rate_vs_paper.csv',rate)
    plt.figure(figsize=(8,4.5))
    for mult in [1,2,3]:
        subset=[r for r in rate if r['transport_area_multiplier_ASSUMED']==mult]
        plt.plot([r['j_mA_cm2'] for r in subset],[r['second_cycle_model_EE_percent'] for r in subset],'-o',label=f'Model, film-area multiplier {mult} (assumed)')
    plt.scatter([30,90],[86.02,64.37],marker='s',s=65,label='Published two exact rate anchors',zorder=5)
    plt.xlabel('Current density (mA/cm²)');plt.ylabel('Energy efficiency (%)')
    plt.title('Model vs published anchors: NOT calibrated, high-current mismatch')
    plt.legend(fontsize=7);fig('13_battery_rate_comparison.png')
    res['rate_sensitivity']=rate
    # Audits
    max_mass=max(abs(r['mass_balance_relative_error']) for r in rfields.values())
    max_cfd=max(abs(r['mass_balance_error']) for r in lobm)
    assert max_mass<1e-8 and max_cfd<1e-8
    res['audit']={'fv_mass_error_max':max_mass,'lbm_coupled_transport_mass_error_max':max_cfd,
                  'porous_felt_FD_verification_error_max':res['porous_felt_independent_FVM_max_error'],
                  'known_mismatch':'3x effective film area gives 90 mA/cmc2 performance but parameter unknown; no independent fit or full pore-scale CFD.'}
    (DATA/'summary.json').write_text(json.dumps(compact(res),indent=2,ensure_ascii=False))
    print('Created',len(list(PLOTS.glob('*.png'))),'figures and',len(list(DATA.glob('*.csv'))),'tables')
    print('Max FV balance error:',res['audit']['fv_mass_error_max'])
    print('LBM hydraulic resistance ratio estimate:',resistance_ratio)
    print('Measured 30/90 EE: 86.02/64.37%. Model results in 11_battery_rate_vs_paper.csv (NOT validated).')

if __name__=='__main__':
    ap=ArgumentParser();ap.add_argument('--fast',action='store_true');args=ap.parse_args()
    main(fast=args.fast)
