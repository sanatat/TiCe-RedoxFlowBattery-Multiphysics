"""Numerical checks, not experimental validation of the real battery."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from transport_fv import Channel, channel_transport,poiseuille_velocity
from electrochem_cell import (Battery, electrolyte_props,init_tanks,metal_inventory,
                              membrane_step,rebalance,faradaic_step,ocv,
                              mass_transport_inputs,voltage,virtual_eis,run_cycles,published_apparent_D)

@pytest.fixture(scope='module')
def km():return mass_transport_inputs(Battery(),nx=64,ny=24)[0]

def test_si_geometry_units():
    b=Battery();assert np.isclose(b.area_m2,9e-4)
    assert np.isclose(b.tank_m3,2e-5)
    assert np.isclose(b.flow_m3_s,20e-6/60)
    assert np.isclose(b.membrane_ASR_ohm_m2,.0942e-4)

def test_si_transport_units_are_apparent():
    assert np.isclose(published_apparent_D('Ti3'),1.78e-8)

def test_fv_zero_wall_reaction_concentration_constant():
    p=Channel()
    r=channel_transport(p,1.78e-8,250.,k_wall_m_s=1e-17,nx=45,ny=20)
    assert np.nanmax(abs(r['c']-250))<1e-6
    assert abs(r['mass_balance_relative_error'])<1e-9

def test_fv_reactive_conserves_moles():
    r=channel_transport(Channel(),1.78e-8,250.,nx=75,ny=27)
    assert abs(r['mass_balance_relative_error'])<1e-10
    assert 0<r['out_mol_s']<r['in_mol_s']
    assert np.isclose(r['in_mol_s']-r['out_mol_s'],r['reacted_mol_s'],rtol=1e-9)
    assert 0<=r['minimum_c_mol_m3']<250

def test_transport_more_diffusion_gives_more_reaction():
    a=channel_transport(Channel(),1.78e-8,250.,nx=70,ny=26)
    b=channel_transport(Channel(),1.78e-10,250.,nx=70,ny=26)
    assert a['reacted_mol_s']>b['reacted_mol_s']

def test_transport_flow_influence():
    a=channel_transport(Channel(flow_m3_s=10e-6/60),1.78e-8,250.,nx=70,ny=26)
    b=channel_transport(Channel(flow_m3_s=40e-6/60),1.78e-8,250.,nx=70,ny=26)
    assert b['reacted_mol_s']>a['reacted_mol_s']
    assert a['in_mol_s']>a['out_mol_s']

def test_fv_mesh_convergence():
    a=channel_transport(Channel(),1.78e-8,250.,nx=62,ny=28)
    b=channel_transport(Channel(),1.78e-8,250.,nx=124,ny=56)
    assert abs(a['k_apparent_m_s']-b['k_apparent_m_s'])/b['k_apparent_m_s']<.12

def test_poiseuille_velocity_integrates_flow():
    p=Channel()
    ux,uy,solid=poiseuille_velocity(p,33,51)
    assert np.isclose((ux[:,0]*(p.height_m/33)*p.width_m).sum(),p.flow_m3_s,rtol=1e-12)
    assert np.all(uy==0)

def test_pump_power_and_viscosity():
    p=Channel();p2=Channel(mu_Pa_s=p.mu_Pa_s*2)
    assert np.isclose(p2.pressure_drop_Pa,2*p.pressure_drop_Pa)
    assert np.isclose(p2.pump_power_two_cells_W,2*p.pump_power_two_cells_W)

def test_electrolyte_concentration_properties():
    s1,m1=electrolyte_props(.5);s2,m2=electrolyte_props(1.0)
    assert np.isclose(s1,24.685)
    assert s1>s2 and m2>m1

def test_nernst_state_of_charge():
    b=Battery()
    assert ocv(init_tanks(b,soc=.8),b)>ocv(init_tanks(b,soc=.2),b)

def test_faradaic_balance():
    b=Battery();c=init_tanks(b)
    before=metal_inventory(c,b)
    fbefore=(c[0,1]+c[1,3])*b.tank_m3
    faradaic_step(c,b,.45,50,True)
    assert np.isclose((c[0,1]+c[1,3])*b.tank_m3-fbefore,2*.45*50/96485.33212)
    assert all(np.isclose(metal_inventory(c,b)[k],v) for k,v in before.items())

def test_membrane_crossover_inventories():
    b=Battery(membrane_ce_m_s=4.54e-9,membrane_ti_m_s=1e-10)
    c=init_tanks(b,soc=.8)
    inv0=metal_inventory(c,b)
    membrane_step(c,b,3600)
    inv=metal_inventory(c,b)
    assert all(abs(inv[k]-v)<1e-14 for k,v in inv0.items())
    assert c[1,1]>0  # some Ce4 moves to negative reservoir

def test_remixing_inventories_and_quench():
    b=Battery();c=init_tanks(b,soc=.8)
    inv0=metal_inventory(c,b)
    rebalance(c)
    assert np.allclose(c[0],c[1])
    assert np.isclose(c[0,1]*c[0,3],0)
    assert all(np.isclose(metal_inventory(c,b)[k],v) for k,v in inv0.items())

def test_electrochemical_voltage_limits(km):
    b=Battery();c=init_tanks(b,soc=.5)
    charge,_,_=voltage(c,b,300,km,True)
    discharge,_,_=voltage(c,b,300,km,False)
    assert np.isfinite(charge) and np.isfinite(discharge)
    assert charge>ocv(c,b)>discharge

def test_eis_physical_impedance(km):
    z=virtual_eis(Battery(),km)
    assert np.all(np.isfinite(z['Z_area_ohm_m2']))
    assert z['Rct_area_ohm_m2']>0 and z['Rdiff_area_ohm_m2']>0
    assert abs(z['Z_area_ohm_m2'][-1].real-z['Rohm_area_ohm_m2'])/z['Rohm_area_ohm_m2']<.2

def test_cycling_mass_conservation(km):
    b=Battery()
    trace,summ,_,inv=run_cycles(b,km,cycles=3,j_mA_cm2=50,regenerate_after=(2,))
    assert len(trace)>10
    assert len(summ)==3
    assert summ[0]['charge_Ah']>0 and summ[0]['discharge_Ah']>0
    assert np.isfinite(summ[1]['EE_percent'])
    for k in ('Ce_mol','Ti_mol'):
        assert abs(inv['initial_inventory_mol'][k]-inv['final_inventory_mol'][k])<1e-11

def test_felt_analytical_numerical_agree():
    from porous_felt import Felt, effectiveness,finite_difference_flux
    f=Felt(active_area_m2_m3=25000.)
    result=finite_difference_flux(f=f,n=480)
    assert result['relative_error']<.005
    assert 0<effectiveness(f=f)['effectiveness']<1

def test_felt_more_electroactive_area_higher_flux():
    from porous_felt import Felt,effectiveness
    a=effectiveness(f=Felt(active_area_m2_m3=2500))
    b=effectiveness(f=Felt(active_area_m2_m3=50000))
    assert b['flux_mol_m2_s']>a['flux_mol_m2_s']
    assert b['effectiveness']<a['effectiveness']

def test_channel_felt_coupling_mass_balance():
    from porous_felt import effectiveness
    kw=effectiveness()['k_wall_m_s']
    r=channel_transport(Channel(),1.78e-8,250,k_wall_m_s=kw,nx=62,ny=28)
    assert abs(r['mass_balance_relative_error'])<1e-10
    assert r['reacted_mol_s']>0

def test_lbm_with_obstacle_conservative_reactive_transport():
    from transport_fv import lbm_velocity
    ch=Channel()
    ux,uy,solid,_=lbm_velocity(ch,nx=64,ny=22,baffles=True,steps=1000)
    a=channel_transport(ch,1.78e-8,250.,flow_field=(ux,uy,solid))
    assert abs(a['mass_balance_relative_error'])<1e-9
    assert a['reacted_mol_s']>0
    assert np.isfinite(a['c'][~solid]).all()

def test_membrane_higher_permeability_more_crossover():
    b1=Battery(membrane_ce_m_s=4.54e-9,membrane_ti_m_s=0.)
    b2=Battery(membrane_ce_m_s=4.54e-8,membrane_ti_m_s=0.)
    c1=init_tanks(b1,soc=.8);c2=init_tanks(b2,soc=.8)
    membrane_step(c1,b1,300);membrane_step(c2,b2,300)
    assert c2[1,1]>c1[1,1]

def test_cycling_lower_current_better_second_cycle_efficiency(km):
    b=Battery()
    _,a,_,_=run_cycles(b,km,cycles=2,j_mA_cm2=30)
    _,c,_,_=run_cycles(b,km,cycles=2,j_mA_cm2=50)
    assert a[1]['EE_percent']>c[1]['EE_percent']

def test_ti_proton_stoichiometry_changes_ideal_formal_voltage():
    from electrolyte_chem import formal_cell_voltage_ideal_proton
    a=formal_cell_voltage_ideal_proton(relative_H_activity=1.)
    bb=formal_cell_voltage_ideal_proton(relative_H_activity=2.)
    expected=-2*8.314462618*298.15/96485.33212*np.log(2)
    assert np.isclose(bb-a,expected)


def test_dft_SI_energy_records_complete():
    from electrolyte_chem import S2
    assert len(S2['Ce3']['binding_Ha'])==9
    assert len(S2['Ce4']['binding_Ha'])==11
    assert S2['Ti4_binuclear']['ligands'][2]==4
