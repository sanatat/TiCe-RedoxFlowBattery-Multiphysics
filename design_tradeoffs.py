"""Conditional flow-field scenario: NOT experimentally validated flow-rate response.
Fitted Fig6c curvature threshold is scaled by the channel-only FV k(Q)/k(20).
Assumes the curvature is entirely transport-controlled, which is NOT identifiable.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from calibrate_transport import discharge_voltage, ROOT

def run():
    out=ROOT/'validation_outputs';figs=ROOT/'validation_figures'
    s=json.loads((out/'validation_summary.json').read_text())
    p=s['fit_parameters'];par=[p['E_extrapolated_V'],p['R_lumped_ohm_cm2'],p['j_curvature_mA_cm2']]
    d=pd.read_csv(out/'fvm_mass_transfer_D_and_Q_sensitivity.csv')
    allrows=[]
    for ds in [1.,.1,.01]:
        sub=d[(d.species=='Ti3')&(d.D_scale==ds)].sort_values('Q_mL_min')
        ref=sub[sub.Q_mL_min==20].iloc[0].k_channel_m_s
        for _,r in sub.iterrows():
            mult=r.k_channel_m_s/ref
            p2=[par[0],par[1],par[2]*mult]
            for jj in [50.,90.,150.,250.]:
                v=discharge_voltage(np.array([jj]),p2)[0]
                allrows.append({'D_scale':ds,'flow_mL_min':r.Q_mL_min,'j_mA_cm2':jj,
                                'assumed_k_flow_ratio':mult,
                                'assumed_threshold_mA_cm2':p2[2],
                                'voltage_conditional_V':float(v) if np.isfinite(v) else np.nan,
                                'power_density_conditional_mW_cm2':float(v*jj) if np.isfinite(v) else np.nan,
                                'channel_only_pump_W_illustrative':r.pump_power_two_halves_W_assumed_geometry,
                                'source_quality':'conditional illustration, no Q-scan experimental validation'})
    result=pd.DataFrame(allrows)
    result.to_csv(out/'conditional_flow_tradeoffs.csv',index=False)
    fig,axs=plt.subplots(1,2,figsize=(12,4.5))
    q=result[(result.D_scale==1)&(result.j_mA_cm2==150)]
    axs[0].plot(q.flow_mL_min,q.voltage_conditional_V,marker='o')
    axs[0].set(xlabel='Flow rate, mL/min',ylabel='Predicted discharge V at 150 mA/cm²',title='Conditional transport-limited flow improvement')
    axs[1].plot(q.flow_mL_min,q.channel_only_pump_W_illustrative*1e3,marker='o')
    axs[1].set(xlabel='Flow rate, mL/min',ylabel='Assumed smooth-channel pumping power, mW',title='Pumping excludes porous-felt loss')
    for a in axs:a.grid(alpha=.25)
    fig.suptitle('NOT validated to reported flow-rate experiments; Fig6c reports one operating flow')
    fig.tight_layout();fig.savefig(figs/'08_conditional_flow_field_tradeoff.png',dpi=160);plt.close(fig)
    return result

if __name__=='__main__':
    out=run();print(out[(out.D_scale==1)&(out.j_mA_cm2==150)].to_string(index=False))
