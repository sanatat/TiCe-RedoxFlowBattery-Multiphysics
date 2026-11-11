"""Conditional cell-stack scaling with parallel hydraulic manifold branches.

The cells are electrically in series and hydraulically in parallel. Equal
pressure across manifold branches and a *linear Darcy* hydraulic law with
assumed branch-resistance spreads yield unequal flows. This is an
engineering SCENARIO, not a tested Ti/Ce module or manufacturer stack model.
Shunt currents, manifolds' spatial losses, thermal gradients, and leakage are
NOT represented.
"""
from dataclasses import dataclass,replace
import numpy as np
from extensions.reactor import ReactorScenario,simulate_discharge


@dataclass(frozen=True)
class StackScenario:
    n_cells:int=8
    current_density_mA_cm2:float=50.0
    nominal_flow_mL_min_per_cell:float=20.0
    hydraulic_resistance_spread:float=0.15  # half span of Darcy branch resistance, hypothetical
    pump_efficiency:float=0.55
    soc:float=0.5
    def __post_init__(self):
        if self.n_cells<1 or min(self.current_density_mA_cm2,self.nominal_flow_mL_min_per_cell,self.pump_efficiency)<=0:
            raise ValueError('invalid stack configuration')
        if not 0<=self.hydraulic_resistance_spread<1:
            raise ValueError('resistance spread must be in [0,1)')


def stack_operating_point(s:StackScenario, base:ReactorScenario=ReactorScenario()):
    n=s.n_cells
    factors=(np.linspace(1-s.hydraulic_resistance_spread,1+s.hydraulic_resistance_spread,n)
             if n>1 else np.ones(1))
    # Darcy hydraulic branches obey DeltaP_i = K_i*Q_i, not Q_i**2.
    # Values below deliberately omit manifold entrance/nonlinear losses.
    weights_lin=1/factors
    flows_lin=s.nominal_flow_mL_min_per_cell*n*weights_lin/weights_lin.sum()
    cells=[simulate_discharge(replace(base,j_mA_cm2=s.current_density_mA_cm2,
                                      flow_mL_min=float(q),soc=s.soc,
                                      pump_efficiency=s.pump_efficiency)) for q in flows_lin]
    if any(not c['feasible'] for c in cells):
        return {'feasible':False,'n_cells':n,'individual_flows_mL_min':flows_lin.tolist(),
                'reason':'Low flow branch exceeds modelled local limiting current'}
    delta_p=np.array([c['pressure_drop_Pa']*factor for c,factor in zip(cells,factors)])
    assert np.allclose(delta_p,delta_p[0],rtol=1e-12)
    q_total_m3_s=s.nominal_flow_mL_min_per_cell*n*1e-6/60
    # Two electrolyte loops; same per-loop total flow and pressure.
    pumping_W=2*delta_p[0]*q_total_m3_s/s.pump_efficiency
    current_A=cells[0]['current_A']
    voltages=np.array([c['cell_voltage_V'] for c in cells])
    power_W=current_A*voltages.sum()
    return {'feasible':True,'n_cells':n,
            'individual_flows_mL_min':flows_lin.tolist(),
            'flow_balance_relative_error':abs(flows_lin.sum()-s.nominal_flow_mL_min_per_cell*n)/
                                           (s.nominal_flow_mL_min_per_cell*n),
            'branch_pressure_spread_Pa':float(np.ptp(delta_p)),
            'stack_voltage_V':float(voltages.sum()),'min_cell_voltage_V':float(voltages.min()),
            'max_cell_voltage_V':float(voltages.max()),'gross_power_W':float(power_W),
            'pump_power_W':float(pumping_W),'net_power_W':float(power_W-pumping_W),
            'pump_fraction_of_gross':float(pumping_W/power_W) if power_W>0 else float('nan'),
            'cell_voltages_V':voltages.tolist(),
            'caveat':'conditional Darcy hydraulic network; neglects shunt currents and cell nonuniformity'}
