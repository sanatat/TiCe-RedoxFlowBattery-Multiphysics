"""Numerical verification, not experimental validation."""
from dataclasses import replace
import math
import numpy as np
import pytest
from extensions.chemistry import MixedAcidHypothesis,conditional_partition,mixed_acid_scenario
from extensions.reactor import ReactorScenario,simulate_discharge
from extensions.membrane import MembraneScenario,simulate_membrane
from extensions.stack import StackScenario,stack_operating_point
from extensions.lifetime import LifetimeScenario,simulate_lifetime


def test_conditional_partition_mass_conservation():
    for beta in (0,1,5,50):
        r=conditional_partition(0.5,1.0,beta)
        assert r['free_M']+r['bound_M']==pytest.approx(.5)
        assert 0<=r['free_M']<=.5

def test_chemistry_msa_anchor():
    r=mixed_acid_scenario(MixedAcidHypothesis(sulfate_fraction=0),.5)
    assert r['conductivity_mS_cm_scenario']==pytest.approx(246.85)
    assert r['viscosity_cP_scenario']==pytest.approx(2.76)
    assert not r['thermodynamically_validated']

def test_acid_fraction_bounds():
    with pytest.raises(ValueError):MixedAcidHypothesis(sulfate_fraction=1.1)

def test_acid_zero_beta_no_complexation():
    r=mixed_acid_scenario(MixedAcidHypothesis(conditional_beta_Ce_M_inv=0),.5)
    assert r['Ce_bound_M_conditional']==pytest.approx(0)

def test_reactor_balances_both_elements():
    r=simulate_discharge(ReactorScenario())
    assert r['feasible']
    assert r['Ce_balance_relative_error']<1e-10
    assert r['Ti_balance_relative_error']<1e-10
    assert r['species_reacted_mol_s']==pytest.approx(r['current_A']/96485.33212)

def test_reactor_voltage_falls_with_current():
    vs=[simulate_discharge(ReactorScenario(j_mA_cm2=j))['cell_voltage_V'] for j in (30,50,90)]
    assert vs[0]>vs[1]>vs[2]

def test_reactor_lower_flow_higher_single_pass_utilization():
    a=simulate_discharge(ReactorScenario(flow_mL_min=20))
    b=simulate_discharge(ReactorScenario(flow_mL_min=30))
    assert a['utilization_single_pass']>b['utilization_single_pass']
    assert a['pressure_drop_Pa']<b['pressure_drop_Pa']

def test_reactor_thicker_felt_reduces_pressure_drop_only_under_darcy():
    a=simulate_discharge(ReactorScenario(felt_thickness_m=.005))
    b=simulate_discharge(ReactorScenario(felt_thickness_m=.010))
    assert b['pressure_drop_Pa']<a['pressure_drop_Pa']

def test_low_effective_D_can_limit_current():
    r=simulate_discharge(ReactorScenario(j_mA_cm2=90,apparent_D_scale=.05))
    assert not r['feasible']

def test_membrane_inventory_and_volume_conservation():
    s=MembraneScenario(water_transfer_m3_s=1e-12,z_Ce_effective=-1,potential_drop_V=.01)
    r=simulate_membrane(s,hours=100)
    assert r['final_Ce_mol']==pytest.approx(r['initial_Ce_mol'],rel=1e-12)
    assert r['final_Ti_mol']==pytest.approx(r['initial_Ti_mol'],rel=1e-12)
    assert r['final_volume_m3']==pytest.approx(r['initial_volume_m3'])

def test_membrane_equal_concentrations_no_field_no_water_stays_equal():
    s=MembraneScenario()
    r=simulate_membrane(s,hours=10,initial_concentrations_mol_m3=((500,200),(500,200)))
    assert r['trace'][-1]['Ce_A_M']==pytest.approx(.5)
    assert r['trace'][-1]['Ti_B_M']==pytest.approx(.2)

def test_membrane_120h_magnitude_is_scenario_not_validation():
    s=MembraneScenario(Ti_permeability_m_s=0)
    r=simulate_membrane(s,hours=120,initial_concentrations_mol_m3=((500,0),(0,0)))
    assert .03 < r['trace'][-1]['Ce_B_M'] <.06

def test_stack_one_cell_consistency():
    r=stack_operating_point(StackScenario(n_cells=1,hydraulic_resistance_spread=0))
    b=simulate_discharge(ReactorScenario())
    assert r['feasible']
    assert r['stack_voltage_V']==pytest.approx(b['cell_voltage_V'])

def test_stack_series_voltage_and_flow_balance():
    r=stack_operating_point(StackScenario(n_cells=6,hydraulic_resistance_spread=0))
    single=stack_operating_point(StackScenario(n_cells=1,hydraulic_resistance_spread=0))
    assert r['stack_voltage_V']==pytest.approx(single['stack_voltage_V']*6)
    assert r['flow_balance_relative_error']<1e-12
    assert r['branch_pressure_spread_Pa']<1e-8

def test_stack_maldistribution():
    r=stack_operating_point(StackScenario(n_cells=6,hydraulic_resistance_spread=.2))
    assert r['feasible']
    assert max(r['individual_flows_mL_min'])>min(r['individual_flows_mL_min'])
    assert r['flow_balance_relative_error']<1e-12

def test_lifetime_no_remix_monotone_loss():
    r=simulate_lifetime(LifetimeScenario(remix_enabled=False))
    c=[v['capacity_fraction'] for v in r['trace']]
    assert all(x>=y for x,y in zip(c,c[1:]))

def test_lifetime_remix_recovers_only_reversible_component():
    r=simulate_lifetime(LifetimeScenario())
    assert len(r['remix_cycles'])>0
    for row in r['trace']:
        if row['remixed']:
            assert row['capacity_fraction']>row['pre_remix_capacity_fraction']
            assert row['capacity_fraction']<=row['felt_activity_fraction']*row['active_inventory_fraction']+1e-12

def test_lifetime_irreversible_deterioration_survives_remix():
    r=simulate_lifetime(LifetimeScenario())
    assert r['final_felt_activity_fraction']<1
    assert r['final_active_inventory_fraction']<1
