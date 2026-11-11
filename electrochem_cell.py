"""Ti-Ce mixed active species flow-battery reservoir-and-cell demonstrator.

Measured SI geometry + reported apparent kinetic/transport values; contact
resistance, fluid-channel dimensions, interfacial Cdl, local side chemistry,
membrane Ti permeability and pump efficiency are illustrative assumptions.
Finite-volume 2D channel wall-transfer coefficients drive cell overpotentials.
This is not an experimentally validated battery digital twin.
"""
from dataclasses import dataclass
from pathlib import Path
import json
import numpy as np
from scipy.integrate import trapezoid
from transport_fv import Channel,channel_transport

F=96485.33212;R=8.314462618
SPECIES=('Ce3','Ce4','Ti4','Ti3')
SI=json.loads((Path(__file__).parent/'si_reference.json').read_text())

@dataclass(frozen=True)
class Battery:
    area_m2:float=SI['cell']['area_cm2']*1e-4
    tank_m3:float=SI['cell']['reservoir_each_mL']*1e-6
    flow_m3_s:float=SI['cell']['flow_each_mL_min']*1e-6/60
    cm_mol_m3:float=500.0
    temp_K:float=298.15
    acid_mol_L:float=2.0
    membrane_ASR_ohm_m2:float=SI['membrane']['membrane_ASR_before_ohm_cm2']*1e-4
    membrane_ce_m_s:float=SI['membrane']['Hcell_Ce_permeability_cm_s']*1e-2
    membrane_ti_m_s:float=2e-10  # unknown, scenario assumption; 'not detected' is not zero
    extra_ASR_ohm_m2:float=2.0e-4 # unknown cell/contact/porous resistances
    bulk_ion_path_m:float=0.0005 # illustrative effective ionic conduction path
    active_area_factor:float=1.0 # effective heterogeneous reaction surface relative to geo
    quench_m3_mol_s:float=0.0 # unknown homogeneous Ce4+ + Ti3 quench kinetics
    max_charge_s:float=1800 # paper reported 30min or upper cutoff
    cutoff_charge_V:float=SI['cell']['charge_cutoff_V']
    cutoff_discharge_V:float=SI['cell']['discharge_cutoff_V']
    pump_efficiency:float=0.55
    @property
    def current_mA_cm2_to_A(self):return self.area_m2*10.0


def electrolyte_props(ce_ti_each_M=.5,T_C=None):
    conc=SI['conductivity_viscosity_vs_concentration']
    a=np.asarray(conc['M_each_metal'],float)
    sigma=np.interp(ce_ti_each_M,a,conc['conductivity_mS_cm'])*.1  # mS/cm -> S/m
    visc=np.interp(ce_ti_each_M,a,conc['viscosity_cP'])*1e-3
    if T_C is not None:
        temp=SI['conductivity_viscosity_vs_temperature']
        # 'room' is not precisely stated; take 25C as illustrative interpolation origin
        tt=np.array([25,30,35,45,65],float)
        sT=np.interp(T_C,tt,temp['conductivity_mS_cm'])*.1
        muT=np.interp(T_C,tt,temp['viscosity_cP'])*1e-3
        sref=temp['conductivity_mS_cm'][0]*.1
        muref=temp['viscosity_cP'][0]*1e-3
        sigma*=sT/sref;visc*=muT/muref
    return float(sigma),float(visc)


def channel_for_battery(b:Battery,concentration_M=.5,T_C=None,height_m=.001):
    _,mu=electrolyte_props(concentration_M,T_C)
    return Channel(length_m=.03,height_m=height_m,width_m=.03,
                   flow_m3_s=b.flow_m3_s,mu_Pa_s=mu,efficiency=b.pump_efficiency)


def published_apparent_D(name):return SI['electrochemical']['D_mixed_cm2_s'][name]*1e-4

def published_apparent_k(name):return SI['electrochemical']['k_mixed_cm_s'][name]*1e-2


def mass_transport_inputs(b:Battery,d_scale=1.0,channel_height_m=.001,nx=88,ny=30):
    """Return 4 convective/diffusive channel wall conductances for battery.
    Effective wall is perfectly absorbing to estimate an upper-bound
    mass-transfer conductance; porous felt roughness is NOT reproduced.
    """
    ch=channel_for_battery(b,height_m=channel_height_m)
    r={s:channel_transport(ch,published_apparent_D(s)*d_scale,
                           250.,k_wall_m_s=1.,nx=nx,ny=ny)
       for s in SPECIES}
    return {s:r[s]['k_apparent_m_s'] for s in SPECIES},r


def init_tanks(b:Battery,soc=.2,mixed=True):
    """2 x 4 reservoirs, all concentrations mol/m3.
    Both reservoirs contain 0.5M Ce and Ti in mixed configuration.
    Spectator ions are placed as inactive starting valence states.
    """
    c=np.zeros((2,4),float)
    c[0,0]=b.cm_mol_m3*(1-soc);c[0,1]=b.cm_mol_m3*soc
    c[1,2]=b.cm_mol_m3*(1-soc);c[1,3]=b.cm_mol_m3*soc
    if mixed:
        c[0,2]=b.cm_mol_m3;c[1,0]=b.cm_mol_m3
    return c


def metal_inventory(c,b:Battery):
    return {'Ce_mol':float(np.sum(c[:,0:2])*b.tank_m3),
            'Ti_mol':float(np.sum(c[:,2:4])*b.tank_m3)}


def ocv(c,b:Battery):
    eps=1e-9
    ce3,ce4=c[0,:2];ti4,ti3=c[1,2:]
    return 1.51+(R*b.temp_K/F)*np.log(max(ce4,eps)/max(ce3,eps)*max(ti3,eps)/max(ti4,eps))


def electrode_potential(E0,c_red,c_ox,j_oxidation_A_m2,k0,T):
    """Exact symmetric Butler-Volmer inversion for solution red -> ox + e-.
    Avoids the singular low-SOC j0 approximation, so charging a chemically
    reset electrolyte is feasible. Electrode potentials relative SHE.
    """
    red=max(0.,float(c_red));ox=max(0.,float(c_ox))
    b=j_oxidation_A_m2/(F*max(k0,1e-20))
    disc=np.hypot(b,2*np.sqrt(red*ox))
    if b>=0: z=(b+disc)/(2*max(red,1e-15))
    else: z=2*max(ox,1e-15)/(disc-b)
    return E0+2*R*T/F*np.log(max(z,1e-90))


def loss_terms(c,b:Battery,j_A_m2:float,km:dict,charge=False):
    """Potential losses under an effective BV + ohmic + film model.
    Reported SI k interpreted as effective heterogeneous rate constants.
    The geometry-specific reactive surface ratio remains unknown.
    """
    ce3,ce4=c[0,:2];ti4,ti3=c[1,2:]
    ce_react=ce3 if charge else ce4
    ti_react=ti4 if charge else ti3
    k_ce=np.sqrt(published_apparent_k('Ce3')*published_apparent_k('Ce4'))*b.active_area_factor
    k_ti=np.sqrt(published_apparent_k('Ti3')*published_apparent_k('Ti4'))*b.active_area_factor
    sign=1. if charge else -1.
    ece=electrode_potential(1.61,ce3,ce4,sign*j_A_m2,k_ce,b.temp_K)
    eti=electrode_potential(.10,ti3,ti4,-sign*j_A_m2,k_ti,b.temp_K)
    e_kin=ece-eti
    E_eq=ocv(c,b)
    eta_activation=sign*(e_kin-E_eq)
    conductivity,_=electrolyte_props(b.cm_mol_m3/1000)
    eta_ohm_mem=j_A_m2*b.membrane_ASR_ohm_m2
    eta_ohm_bulk=j_A_m2*b.bulk_ion_path_m/conductivity
    eta_ohm_extra=j_A_m2*b.extra_ASR_ohm_m2
    jlim_ce=F*km['Ce3' if charge else 'Ce4']*max(ce_react,0)
    jlim_ti=F*km['Ti4' if charge else 'Ti3']*max(ti_react,0)
    limit=min(jlim_ce,jlim_ti)
    if j_A_m2>=.995*limit:return None,limit
    eta_mt=-(R*b.temp_K/F)*(np.log1p(-j_A_m2/jlim_ce)+np.log1p(-j_A_m2/jlim_ti))
    return {'activation_Ce_Ti_V':float(eta_activation),
            'membrane_ohmic_V':float(eta_ohm_mem),'bulk_electrolyte_V':float(eta_ohm_bulk),
            'extra_contact_V':float(eta_ohm_extra),'concentration_V':float(eta_mt),
            'kinetic_cell_voltage_V':float(e_kin)},limit


def voltage(c,b:Battery,j_A_m2,km,charge):
    terms,limit=loss_terms(c,b,j_A_m2,km,charge)
    if terms is None:return np.nan,None,limit
    sign=1 if charge else -1
    ohm_conc=terms['membrane_ohmic_V']+terms['bulk_electrolyte_V']+terms['extra_contact_V']+terms['concentration_V']
    return terms['kinetic_cell_voltage_V']+sign*ohm_conc,terms,limit


def membrane_step(c,b:Battery,dt):
    """Symmetric Fick transport + optional chemical quench; volume fixed.
    No proton migration, convective water transfer or membrane electrodynamics.
    Fick flux conserves Ce and Ti exactly across two reservoirs.
    """
    p=np.array([b.membrane_ce_m_s]*2+[b.membrane_ti_m_s]*2)
    dc=(c[0]-c[1])*p*b.area_m2*dt/b.tank_m3
    c[0]-=dc;c[1]+=dc
    if b.quench_m3_mol_s>0:
        for side in range(2):
            reacted=min(c[side,1],c[side,3],
                        b.quench_m3_mol_s*c[side,1]*c[side,3]*dt)
            c[side,1]-=reacted;c[side,3]-=reacted
            c[side,0]+=reacted;c[side,2]+=reacted
    return c


def faradaic_step(c,b:Battery,current_A,dt,charge):
    dc=current_A*dt/(F*b.tank_m3)
    s=1 if charge else -1
    c[0,0]-=s*dc;c[0,1]+=s*dc
    c[1,2]-=s*dc;c[1,3]+=s*dc
    return c


def rebalance(c):
    """Idealized mixing of equal tank volumes followed by stoichiometric
    Ce4+ + Ti3 -> Ce3+ + Ti4. Assumes spontaneous reaction completion;
    *not* a measured regeneration kinetic law. Rebalancing resets SOC.
    """
    common=(c[0]+c[1])*.5
    reduced=min(common[1],common[3])
    common[1]-=reduced;common[3]-=reduced
    common[0]+=reduced;common[2]+=reduced
    c[:]=common
    return c


def run_cycles(b:Battery,km,cycles=4,j_mA_cm2=50,dt_s=6.,regenerate_after=(),mixed=True):
    c=init_tanks(b,soc=0.001,mixed=mixed)
    start_inv=metal_inventory(c,b)
    ch=channel_for_battery(b)
    pump_power=ch.pump_power_two_cells_W
    j=j_mA_cm2*10.;I=j*b.area_m2
    rows=[];summary=[];time_s=0
    for cy in range(1,cycles+1):
        records=[]
        for charge in (True,False):
            charge_Ah=0;electric_Wh=0;phase_s=0;status='max_steps'
            for step in range(12000):
                v,terms,jlim=voltage(c,b,j,km,charge)
                if not np.isfinite(v):status='transport_limited';break
                if charge and v>=b.cutoff_charge_V:status='voltage_cutoff';break
                if not charge and v<=b.cutoff_discharge_V:status='voltage_cutoff';break
                if charge and phase_s>=b.max_charge_s:status='charge_timer';break
                limited=min(c[0,0],c[1,2]) if charge else min(c[0,1],c[1,3])
                actual_dt=min(dt_s,.75*limited*b.tank_m3*F/I)
                if charge:actual_dt=min(actual_dt,b.max_charge_s-phase_s)
                if actual_dt<=1e-7:status='reactant_depleted';break
                faradaic_step(c,b,I,actual_dt,charge)
                membrane_step(c,b,actual_dt)
                if c.min() < -1e-6:raise ArithmeticError('negative concentration')
                phase_s+=actual_dt;time_s+=actual_dt
                charge_Ah+=I*actual_dt/3600
                electric_Wh+=v*I*actual_dt/3600
                if step%8==0:
                    rows.append({'time_h':time_s/3600,'cycle':cy,
                                 'phase':'charge' if charge else 'discharge',
                                 'voltage_V':float(v),'passed_Ah':float(charge_Ah),
                                 'Ce4_positive_M':c[0,1]/1000,
                                 'Ti3_negative_M':c[1,3]/1000,
                                 'eta_activation_V':terms['activation_Ce_Ti_V'],
                                 'eta_ohmic_V':terms['membrane_ohmic_V']+terms['bulk_electrolyte_V']+terms['extra_contact_V'],
                                 'eta_concentration_V':terms['concentration_V'],
                                 'limiting_mA_cm2':jlim/10})
            records.append((charge_Ah,electric_Wh,phase_s,status))
        chg,dis=records
        ce=dis[0]/chg[0] if chg[0]>0 else np.nan
        ee=dis[1]/chg[1] if chg[1]>0 else np.nan
        pump_ch=pump_power*chg[2]/3600; pump_dis=pump_power*dis[2]/3600
        ee_with_pump=(dis[1]-pump_dis)/(chg[1]+pump_ch) if chg[1]>0 else np.nan
        summary.append({'cycle':cy,'charge_Ah':chg[0],'discharge_Ah':dis[0],
                        'CE_percent':100*ce,'VE_percent':100*ee/ce if ce>0 else np.nan,
                        'EE_percent':100*ee,'EE_with_pump_percent':100*ee_with_pump,
                        'charge_end':chg[3],'discharge_end':dis[3],
                        'regenerated_after':cy in regenerate_after,'time_h':time_s/3600})
        if cy in regenerate_after: rebalance(c)
    final_inv=metal_inventory(c,b)
    return rows,summary,c,{'initial_inventory_mol':start_inv,'final_inventory_mol':final_inv,
                            'pump_power_two_cells_W':pump_power}


def virtual_eis(b:Battery,km,soc=.5,cdl_F_m2=50.,n=140):
    """Small-signal equivalent circuit INFERRED from 0D loss derivatives.
    Synthetic, not measured impedance. Includes assumed capacitance.
    """
    c=init_tanks(b,soc=soc)
    base,_=loss_terms(c,b,1e-5,km,False)
    rct=F**-1*R*b.temp_K*(1/(F*np.sqrt(published_apparent_k('Ce3')*published_apparent_k('Ce4'))*np.sqrt(c[0,0]*c[0,1])*b.active_area_factor)
                               +1/(F*np.sqrt(published_apparent_k('Ti3')*published_apparent_k('Ti4'))*np.sqrt(c[1,2]*c[1,3])*b.active_area_factor))
    sigma,_=electrolyte_props()
    r0=b.membrane_ASR_ohm_m2+b.bulk_ion_path_m/sigma+b.extra_ASR_ohm_m2
    jlim_ce=F*km['Ce4']*c[0,1];jlim_ti=F*km['Ti3']*c[1,3]
    rw=R*b.temp_K/F*(1/jlim_ce+1/jlim_ti)
    delta=published_apparent_D('Ti3')/km['Ti3']
    tau=delta**2/published_apparent_D('Ti3')
    w=2*np.pi*np.logspace(-3,5,n)
    zeta=np.sqrt(1j*w*tau)
    zw=rw*np.tanh(zeta)/zeta
    z=r0+rct/(1+1j*w*rct*cdl_F_m2)+zw
    return {'frequency_Hz':w/(2*np.pi),'Z_area_ohm_m2':z,'Rohm_area_ohm_m2':r0,
            'Rct_area_ohm_m2':rct,'Rdiff_area_ohm_m2':rw,'Cdl_F_m2':cdl_F_m2,
            'diffusion_time_s':tau}
