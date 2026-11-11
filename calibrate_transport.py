"""SI-backed transport inverse study on Fu et al. 2026 Ti-Ce mRFB.

Fig.6c polarization points are *visually digitized approximate readings*
from the published raster image, not original raw experimental data.
The reduced-order fit has one effective, non-identifiable transport threshold;
it MUST NOT be interpreted as a measurement of carbon felt active area.
"""
from pathlib import Path
from dataclasses import replace
import json
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from electrochem_cell import Battery, published_apparent_D, electrolyte_props, channel_for_battery
from transport_fv import channel_transport

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'validation_outputs'
FIG=ROOT/'validation_figures'
FARADAY=96485.33212
R=8.314462618
TEMP=298.15
MEMBRANE_ASR_OHM_CM2=0.0942 # SI Figure S20, before cycling
FIG6C_SIGMA_V=0.04 # estimated visual digitization & line-identification uncertainty, NOT measured SD

# Approximate manual readings from Figure 6(c), blue/cyan discharge 0.50M mRFB branch.
# Units j:mA/cm2, V:V. Underlying Figure 6c raster has crossing curves and low resolution.
# These should NOT be mistaken for primary numerical data supplied by authors.
FIG6C_J=[15,35,55,75,95,115,135,155,175,195,215,235,255,275]
FIG6C_V=[1.68,1.58,1.48,1.38,1.28,1.18,1.09,1.00,0.92,0.85,0.78,0.70,0.56,0.30]
HOLDOUT=[1,4,7,10]


def discharge_voltage(j, p, temp_K=TEMP):
    """Lumped discharge model (illustrative, not mechanistic component separation).

    V = E* - j(ASR_N115+R_lumped) + (2RT/F) ln(1-j/j_lim)
    j in mA/cm2, ASR in Ohm cm2. R_lumped also absorbs linearized activation,
    bulk/contacts and unmodelled losses. j_lim is an effective curvature scale;
    the branch endpoint need NOT equal a true diffusion-limited current.
    """
    j=np.asarray(j,float)
    E0,R_lumped,j_lim=p
    if np.any(j>=j_lim):
        return np.full_like(j,np.nan)
    thermal=2*R*temp_K/FARADAY
    return E0 - j*.001*(R_lumped+MEMBRANE_ASR_OHM_CM2) + thermal*np.log1p(-j/j_lim)


def fit_reduced(j,V,mask):
    lower=[1.48,0.,max(j)+.5]
    upper=[1.95,12.,1e5]
    res=least_squares(lambda p: discharge_voltage(j[mask],p)-V[mask],
                      x0=[1.72,4.,max(j)+20],bounds=(lower,upper),
                      max_nfev=1000,xtol=1e-11,ftol=1e-11,gtol=1e-11)
    return res.x


def load_digitized():
    return pd.DataFrame({'j_mA_cm2':FIG6C_J,'V_discharge_approx_V':FIG6C_V,
                         'approx_digitization_uncertainty_V':FIG6C_SIGMA_V,
                         'held_out_for_assessment':[i in HOLDOUT for i in range(len(FIG6C_J))],
                         'source':['Fu et al 2026 main Fig.6c, blue 0.50M mRFB discharge; visual estimate']*len(FIG6C_J)})


def transport_sensitivity():
    """Actual 2D finite-volume flow/transport solves, assumed 1mm channel geometry.

    Reported Table S6 'D' inferred from porous-felt quasi-reversible CV is apparent.
    D factors are scenarios, not corrections of reported values.
    """
    b=Battery()
    data=[]
    flows=(10,20,40)
    dscales=(0.01,0.1,1.)
    for species in ('Ti3','Ce4'):
        for dscale in dscales:
            for flow in flows:
                ch=replace(channel_for_battery(b),flow_m3_s=flow*1e-6/60)
                res=channel_transport(ch,published_apparent_D(species)*dscale,250.,
                                      k_wall_m_s=1.,nx=72,ny=28)
                data.append({'species':species,'D_scale':dscale,'Q_mL_min':flow,
                             'k_channel_m_s':res['k_apparent_m_s'],
                             'jlim_channel_at_50pctSOC_mA_cm2':FARADAY*res['k_apparent_m_s']*250./10.,
                             'concentration_outlet_min_mol_m3':res['minimum_c_mol_m3'],
                             'mass_balance_relative_error':res['mass_balance_relative_error'],
                             'delta_p_channel_Pa_assumed_geometry':res['delta_p_Pa'],
                             'pump_power_two_halves_W_assumed_geometry':res['pump_power_two_cells_W']})
    return pd.DataFrame(data)


def effective_transport_envelope(p,df_transport):
    """Required effective k for the fitted curvature scale at specified SOC.

    Includes combined uncertainties in actual SOC, potential kinetics,
    channel topology and porous area. NOT a calibrated felt active surface.
    """
    rows=[]
    jlim=p[2]
    for soc in [0.25,0.50,0.75,1.00]:
        c=500.*soc # active species in the electrode at unknown polarization SOC
        k_req=(jlim*10.)/(FARADAY*c)
        for dscale in [0.01,0.1,1.0]:
            for species in ['Ti3','Ce4']:
                kt=df_transport[(df_transport.species==species)&
                                (df_transport.D_scale==dscale)&
                                (df_transport.Q_mL_min==20)].iloc[0].k_channel_m_s
                rows.append({'species':species,'SOC_assumed':soc,'D_scale':dscale,
                             'jlim_fitted_effective_mA_cm2':jlim,
                             'k_required_m_s':k_req,'channel_k_proxy_m_s':kt,
                             'required_effective_multiplier_over_assumed_channel':k_req/kt,
                             'interpretation':'composite wetting/flow topology/interfacial area/kinetics, NOT intrinsic felt area'})
    return pd.DataFrame(rows)


def felt_structural_sensitivity():
    """Illustrative active surface density sensitivities. No unique inversion.

    SI apparent k is NOT used as true intrinsic surface kinetic constant;
    the illustrative k_surface is a separate unknown sensitivity variable.
    """
    data=[]
    for ds in [0.01,.1,1.]:
        for a in [500,2000,10000,25000,100000,500000]:
            for k_surface in [1e-7,1e-6,1e-5]:
                # Effectiveness() currently uses article's apparent k,
                # so here compute a truly illustrative, explicitly separate k.
                D=published_apparent_D('Ti3')*ds*0.85/2.5
                lam=a*k_surface
                L=.005
                phi=L*np.sqrt(lam/D)
                eta=np.tanh(phi)/phi
                kw=np.sqrt(D*lam)*np.tanh(phi)
                data.append({'D_apparent_multiplier':ds,'a_s_m2_m3_assumed':a,
                             'intrinsic_k_surface_m_s_assumed':k_surface,
                             'Thiele':phi,'eta_felt':eta,'k_felt_proxy_m_s':kw,
                             'jlim_proxy_at_c250_mA_cm2':FARADAY*kw*250/10})
    return pd.DataFrame(data)


def bootstrap_digitization(j,v,train,n=240,seed=42):
    rng=np.random.default_rng(seed)
    ps=[]
    for _ in range(n):
        vi=v+rng.normal(0,FIG6C_SIGMA_V,size=len(v))
        ps.append(fit_reduced(j,vi,train))
    return np.array(ps)


def create_plots(df,p,boots,transport,envelope,felt):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    j=df.j_mA_cm2.to_numpy(); v=df.V_discharge_approx_V.to_numpy()
    train=~df.held_out_for_assessment.to_numpy()
    x=np.linspace(0,275,380)
    # 1 digitized + fitted with honest uncertainty
    fig,ax=plt.subplots(figsize=(8,5))
    ax.errorbar(j[train],v[train],yerr=.04,fmt='o',ms=4,capsize=2,label='Approx. Fig6c train (±0.04 V)')
    ax.errorbar(j[~train],v[~train],yerr=.04,fmt='s',ms=5,capsize=2,label='Approx. Fig6c withheld')
    ax.plot(x,discharge_voltage(x,p),label='3-param lumped fit',lw=2)
    ax.set(xlabel='Current density (mA/cm²)',ylabel='Cell discharge voltage (V)',title='Fu et al. Fig.6c – visual readout, NOT raw experimental data')
    ax.legend(fontsize=8);ax.grid(alpha=.2);fig.tight_layout();fig.savefig(FIG/'01_fig6c_constrained_fit.png',dpi=170);plt.close(fig)
    # bootstrap fit p ranges: can demonstrate j_lim non-identifiability
    fig,axs=plt.subplots(1,3,figsize=(12,3.5))
    for k,(lab,unit) in enumerate([('Extrapolated intercept','V'),('Lumped polarization ASR','Ω·cm²'),('Effective curvature threshold','mA/cm²')]):
        axs[k].hist(boots[:,k],bins=25)
        axs[k].axvline(p[k],ls='--',c='black',lw=1)
        axs[k].set(xlabel=f'{lab} ({unit})',ylabel='Simulations')
    fig.suptitle('Bootstrap under assumed ±0.04 V manual-reading uncertainty')
    fig.tight_layout();fig.savefig(FIG/'02_digitization_uncertainty.png',dpi=170);plt.close(fig)
    # transport D sensitivity at 20mL/min
    fig,ax=plt.subplots(figsize=(7,4.5))
    for sp in ['Ti3','Ce4']:
        for dscale,sty in [(1,'-'),(.1,'--'),(.01,':')]:
            sub=transport[(transport.species==sp)&(transport.D_scale==dscale)].sort_values('Q_mL_min')
            ax.plot(sub.Q_mL_min,sub.jlim_channel_at_50pctSOC_mA_cm2,sty,marker='o',label=f'{sp}, D×{dscale:g}')
    ax.set(xlabel='Flow each half-cell (mL/min)',ylabel='Channel-only limiting-current proxy (mA/cm²)',title='2D finite-volume advection–diffusion; ASSUMED 1-mm channel')
    ax.legend(fontsize=8,ncol=2);ax.grid(alpha=.2);fig.tight_layout();fig.savefig(FIG/'03_flow_diffusion_sensitivity.png',dpi=170);plt.close(fig)
    # multipliers for SOC
    fig,ax=plt.subplots(figsize=(8,4.5))
    for dscale in [1,.1,.01]:
        sub=envelope[(envelope.species=='Ti3')&(envelope.D_scale==dscale)]
        ax.plot(sub.SOC_assumed,sub.required_effective_multiplier_over_assumed_channel,marker='o',label=f'Apparent D×{dscale:g}')
    ax.axhline(1,color='gray',ls=':')
    ax.set_yscale('log');ax.set(xlabel='Assumed limiting electrode SOC / available species fraction',ylabel='Inferred k_eff / assumed 1-mm channel k',title='Composite enhancement required for fitted j*; NOT a_s')
    ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(FIG/'04_identifiability_vs_SOC.png',dpi=170);plt.close(fig)
    # porous sensitivity
    fig,ax=plt.subplots(figsize=(8,5))
    for ds in [.01,.1,1.]:
        sub=felt[(felt.D_apparent_multiplier==ds)&(felt.intrinsic_k_surface_m_s_assumed==1e-6)]
        ax.plot(sub.a_s_m2_m3_assumed,sub.jlim_proxy_at_c250_mA_cm2,marker='o',label=f'D apparent ×{ds:g}')
    ax.set_xscale('log');ax.set_yscale('log')
    ax.set(xlabel='Assumed felt specific active area (m²/m³)',ylabel='Illustrative felt-only capacity (mA/cm²)',title='Porous-felt Thiele model: illustrative k_surface=10⁻⁶ m/s')
    ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(FIG/'05_felt_area_families.png',dpi=170);plt.close(fig)
    # losses ASR membrane vs lumped
    fig,ax=plt.subplots(figsize=(7,4))
    currents=np.linspace(0,120,100)
    ax.plot(currents,currents*.001*MEMBRANE_ASR_OHM_CM2,label='Measured N115 ASR term')
    ax.plot(currents,currents*.001*p[1],label='Fitted lumped losses (not pure ohmic)')
    ax.set(xlabel='Current density (mA/cm²)',ylabel='Voltage-loss component (V)',title='The membrane ASR alone cannot explain full polarization')
    ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(FIG/'06_membrane_vs_lumped.png',dpi=170);plt.close(fig)


def run():
    OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
    df=load_digitized();df.to_csv(OUT/'digitized_fig6c_approx.csv',index=False)
    j=df.j_mA_cm2.to_numpy();v=df.V_discharge_approx_V.to_numpy()
    train=~df.held_out_for_assessment.to_numpy()
    p=fit_reduced(j,v,train)
    predict=discharge_voltage(j,p)
    train_rmse=float(np.sqrt(np.mean((predict[train]-v[train])**2)))
    holdout_rmse=float(np.sqrt(np.mean((predict[~train]-v[~train])**2)))
    df['model_voltage_V']=predict;df['residual_V']=predict-v
    df.to_csv(OUT/'fig6c_with_fit_and_holdouts.csv',index=False)
    # Identifiability stress-test: remove the only manually read high-current end point.
    p_without_end=fit_reduced(j,v,train&(j<265))
    end_extrap=float(discharge_voltage(np.array([275.]),p_without_end)[0])
    boots=bootstrap_digitization(j,v,train)
    boot_summary={key:dict(zip(['p05','median','p95'],np.percentile(boots[:,idx],[5,50,95]).tolist()))
                  for idx,key in enumerate(['E_extrapolated_V','R_lumped_ohm_cm2','j_curvature_mA_cm2'])}
    pd.DataFrame(boots,columns=['E_extrapolated_V','R_lumped_ohm_cm2','j_curvature_mA_cm2']).to_csv(OUT/'bootstrap_parameters.csv',index=False)
    transport=transport_sensitivity();transport.to_csv(OUT/'fvm_mass_transfer_D_and_Q_sensitivity.csv',index=False)
    env=effective_transport_envelope(p,transport);env.to_csv(OUT/'nonidentifiable_transport_envelope.csv',index=False)
    felt=felt_structural_sensitivity();felt.to_csv(OUT/'felt_surface_parameter_families.csv',index=False)
    create_plots(df,p,boots,transport,env,felt)
    # High-current endpoint sensitivity and lack of unique j_lim estimation.
    fig,ax=plt.subplots(figsize=(7.8,4.5))
    xx=np.linspace(10,275,350)
    ax.plot(j,v,'o',color='gray',label='Visual Fig.6c estimates')
    ax.plot(xx,discharge_voltage(xx,p),label='Fit with endpoint')
    ax.plot(xx,discharge_voltage(xx,p_without_end),'--',label='Fit without last endpoint')
    ax.set(xlabel='Current density (mA/cm²)',ylabel='Discharge voltage (V)',title='Mass-transfer curvature not independently identifiable')
    ax.grid(alpha=.2);ax.legend();fig.tight_layout();fig.savefig(FIG/'07_endpoint_identifiability.png',dpi=170);plt.close(fig)
    # physically grounded single pass utilization from reported Q and geometric area
    pass_util={}
    for jj in [30,50,90,150,275]:
        current_A=jj*1e-3*9.
        feed_mol_s=(20e-6/60)*500.
        pass_util[str(jj)]=float(current_A/FARADAY/feed_mol_s)
    baseline=transport[(transport.species=='Ti3')&(transport.D_scale==1)&(transport.Q_mL_min==20)].iloc[0]
    summary={'source':'Fu et al. 2026 ChemSusChem 10.1002/cssc.202502732, main Fig6c visually approximated; SI S6-S8/S20',
             'scope':'REDUCED ORDER DISCHARGE ONLY; not full mechanistic multi-parameter calibration',
             'digitization_uncertainty_assumed_V':FIG6C_SIGMA_V,
             'fit_parameters':{'E_extrapolated_V':float(p[0]),'R_lumped_ohm_cm2':float(p[1]),'j_curvature_mA_cm2':float(p[2])},
             'train_rmse_V':train_rmse,'heldout_rmse_V':holdout_rmse,
             'endpoint_removal_stress_test':{'parameters_without_275_point':dict(zip(['E0_V','R_lumped_ohm_cm2','j_curvature_mA_cm2'],map(float,p_without_end))),
                                             'predicted_voltage_at_275_without_point_V':end_extrap,
                                             'manual_275_voltage_V':float(v[-1])},
             'bootstrap_90pct':boot_summary,
             'n_approx_figure_points':len(j),'n_train':int(train.sum()),'n_heldout':int((~train).sum()),
             'j_mA_cm2_visual_curve_endpoint_approx':275.,
             'FVM_assumed_channel_Ti3_k_m_s':float(baseline.k_channel_m_s),
             'FVM_max_conservation_error_abs':float(transport.mass_balance_relative_error.abs().max()),
             'single_pass_utilization_fraction_at_0_5M_and_20mL_min':pass_util,
             'membrane_resistance_ohm_cm2_S20':MEMBRANE_ASR_OHM_CM2,
             'notes':['Membrane ASR is directly measured (S20); lumped R is NOT separate kinetic/ohmic fit',
                      'Only 0.50 M mRFB Fig6c curve fitted, not 0.75 or 1M; no out-of-study validation',
                      'Fitted j_curvature is NOT demonstrated molecular-diffusion-limited current; removing the visually estimated last point makes this parameter effectively unbounded',
                      'Bootstrap apparently narrow j_curvature interval is artificial because fit is constrained above observed max; do not interpret as uncertainty interval for a physical limiting current',
                      'Porous felt active area not identified; wide families of a_s, k_surface and tortuosity fit a given effective transport',
                      'CV-derived D values are unusually high and likely apparent; use 1x,0.1x,0.01x sensitivity',
                      'Missing raw V-I, hydrodynamic geometry, felt surface/porosity, SOC, electrode EIS; validation pending']}
    (OUT/'validation_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf8')
    return summary

if __name__=='__main__':
    s=run()
    print(json.dumps({k:v for k,v in s.items() if k not in ['notes','bootstrap_90pct']},indent=2,ensure_ascii=False))
    print('Wrote',str(OUT),'and',str(FIG))
