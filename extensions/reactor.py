"""Axial channel + porous-electrode reduced-order electrochemical reactor.

At each segment the imposed cell current consumes one electron-equivalent of
Ce(IV) and Ti(III), updates both redox inventories and recomputes local Nernst,
Butler-Volmer, ohmic and film-concentration overpotentials.  It therefore
couples flow, residence/contact, reaction, mass transfer and membrane ASR.

This is NOT a complete 2D/3D Nernst-Planck-Poisson/porous electrode CFD model.
Unknown exchange-current utilization, effective film thickness, felt-specific
area, felt permeability, and contact ASR remain explicit SCENARIO parameters.
"""
from dataclasses import dataclass, replace
import math
import numpy as np
from electrochem_cell import F,R,Battery,electrolyte_props,published_apparent_D,published_apparent_k


@dataclass(frozen=True)
class ReactorScenario:
    j_mA_cm2:float=50.0
    soc:float=0.5
    flow_mL_min:float=20.0
    temp_K:float=298.15
    n_segments:int=60
    felt_thickness_m:float=0.005
    felt_permeability_m2:float=1.e-10  # unmeasured SCENARIO
    surface_utilization_factor:float=0.30 # unmeasured: kCV -> geometric BV mapping
    effective_area_multiplier:float=1.0 # unmeasured porous active area
    film_thickness_m:float=1.5e-4  # unmeasured effective boundary layer
    apparent_D_scale:float=1.0 # SI CV-derived D may not equal bulk molecular D
    extra_ASR_ohm_m2:float=2e-4 # unmeasured other electrical resistance
    pump_efficiency:float=0.55

    def __post_init__(self):
        if not 0<self.soc<1 or self.n_segments<2:raise ValueError('0<SOC<1 and >=2 segments')
        values=(self.j_mA_cm2,self.flow_mL_min,self.temp_K,self.felt_thickness_m,
                self.felt_permeability_m2,self.surface_utilization_factor,
                self.effective_area_multiplier,self.film_thickness_m,
                self.apparent_D_scale,self.pump_efficiency)
        if min(values)<=0 or self.extra_ASR_ohm_m2<0:raise ValueError('positive scenario inputs')


def simulate_discharge(s:ReactorScenario):
    """Steady single-pass discharge at a FIXED SOC at reservoir inlet.

    The two half-cells carry an identical current in series. Local current
    density remains uniform; axial concentrations and overpotentials vary.
    This is a constrained galvanostatic approximation, not a local-current
    redistribution solver. Results are conditional on the scenario inputs.
    """
    area=.0009;length=.03;width=.03
    q=s.flow_mL_min*1e-6/60.0
    I=s.j_mA_cm2*10*area
    sigma,visc=electrolyte_props(.5,s.temp_K-273.15)
    membrane_ASR=Battery().membrane_ASR_ohm_m2
    bulk_ASR=.0005/sigma
    asr=membrane_ASR+bulk_ASR+s.extra_ASR_ohm_m2
    K=s.felt_permeability_m2
    superficial_u=q/(width*s.felt_thickness_m)
    delta_p=visc*length*superficial_u/K
    pump_W=2*delta_p*q/s.pump_efficiency
    ctotal=500.0
    ce4=ctotal*s.soc;ce3=ctotal-ce4
    ti3=ctotal*s.soc;ti4=ctotal-ti3
    xs=[];v_list=[];eta_act=[];eta_mt=[];eta_ohm=[];cce=[];cti=[];limits=[]
    j=I/area
    ti_D=published_apparent_D('Ti3')*s.apparent_D_scale
    ce_D=published_apparent_D('Ce4')*s.apparent_D_scale
    km_ti=ti_D/s.film_thickness_m*s.effective_area_multiplier
    km_ce=ce_D/s.film_thickness_m*s.effective_area_multiplier
    k_ce=math.sqrt(published_apparent_k('Ce3')*published_apparent_k('Ce4'))
    k_ti=math.sqrt(published_apparent_k('Ti3')*published_apparent_k('Ti4'))
    delta_c=I/(F*q*s.n_segments)
    for seg in range(s.n_segments):
        cce.append(ce4);cti.append(ti3)
        jlim_ce=F*km_ce*max(ce4,0)
        jlim_ti=F*km_ti*max(ti3,0)
        limit=min(jlim_ce,jlim_ti)
        limits.append(limit)
        if min(ce4,ti3)<=delta_c or j>=limit:
            # Do not return misleading voltages after mass-transfer depletion
            return {'feasible':False,'limiting_segment':seg,'limiting_j_mA_cm2':limit/10,
                    'scenario':s.__dict__,'pump_power_W':pump_W}
        j0_ce=F*k_ce*math.sqrt(ce3*ce4)*s.surface_utilization_factor*s.effective_area_multiplier
        j0_ti=F*k_ti*math.sqrt(ti3*ti4)*s.surface_utilization_factor*s.effective_area_multiplier
        Eeq=1.51+(R*s.temp_K/F)*math.log((ce4/ce3)*(ti3/ti4))
        activation=2*R*s.temp_K/F*(math.asinh(j/(2*j0_ce))+math.asinh(j/(2*j0_ti)))
        concentration=-R*s.temp_K/F*(math.log1p(-j/jlim_ce)+math.log1p(-j/jlim_ti))
        ohmic=j*asr
        v=Eeq-activation-concentration-ohmic
        xs.append((seg+.5)/s.n_segments*length)
        v_list.append(v);eta_act.append(activation);eta_mt.append(concentration);eta_ohm.append(ohmic)
        ce4-=delta_c;ce3+=delta_c
        ti3-=delta_c;ti4+=delta_c
    out_v=float(np.mean(v_list))
    mass_in_I_equiv=q*(ctotal*s.soc)
    actual_ce=q*(ctotal*s.soc-ce4)
    actual_ti=q*(ctotal*s.soc-ti3)
    return {'feasible':True,'j_mA_cm2':s.j_mA_cm2,'cell_voltage_V':out_v,
            'power_W':I*out_v,'net_power_W':I*out_v-pump_W,
            'current_A':I,'pump_power_W':pump_W,'pressure_drop_Pa':delta_p,
            'geometric_felt_volume_over_flow_s':area*s.felt_thickness_m/q,  # NOT porosity-adjusted residence time
            'species_reacted_mol_s':I/F,'Ce_reacted_mol_s':actual_ce,
            'Ti_reacted_mol_s':actual_ti,'Ce_balance_relative_error':abs(actual_ce-I/F)/max(I/F,1e-30),
            'Ti_balance_relative_error':abs(actual_ti-I/F)/max(I/F,1e-30),
            'utilization_single_pass':(I/F)/mass_in_I_equiv,
            'exit_Ce4_mol_m3':ce4,'exit_Ti3_mol_m3':ti3,
            'mean_activation_V':float(np.mean(eta_act)),
            'mean_concentration_V':float(np.mean(eta_mt)),
            'mean_ohmic_V':float(np.mean(eta_ohm)),
            'minimum_local_limiting_j_mA_cm2':min(limits)/10,
            'x_m':np.array(xs),'voltage_profile_V':np.array(v_list),
            'Ce4_profile_mol_m3':np.array(cce),'Ti3_profile_mol_m3':np.array(cti),
            'scenario':s.__dict__}
