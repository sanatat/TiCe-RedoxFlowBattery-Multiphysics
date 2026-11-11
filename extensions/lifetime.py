"""Scenario-based capacity fade: reversible tank imbalance + irreversible felt fade.

No longitudinal experimental Ti/Ce composition or measured kinetic ageing
constants are available for calibrating these mechanisms. Values are NOT
reproduction of Fu et al.'s 3,721 measured cycles. Mixing restores a *modeled*
reversible component but does not regenerate electrode activity or lost metal.
"""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class LifetimeScenario:
    cycles:int=1000
    imbalance_growth_per_cycle:float=0.00055  # hypothetical reversible loss
    felt_degradation_per_cycle:float=0.00006  # hypothetical irreversible activity fade
    irreversible_inventory_loss_per_cycle:float=0.00001  # hypothetical side-chemistry effect
    remix_threshold_fraction:float=0.80
    recovered_imbalance_fraction:float=0.96  # fraction of imbalance removed at remix
    remix_enabled:bool=True
    def __post_init__(self):
        if self.cycles<1:raise ValueError('cycles>0')
        for x in (self.imbalance_growth_per_cycle,self.felt_degradation_per_cycle,
                  self.irreversible_inventory_loss_per_cycle,self.remix_threshold_fraction,
                  self.recovered_imbalance_fraction):
            if not 0<=x<=1:raise ValueError('rates/fractions in [0,1]')
        if self.remix_threshold_fraction<=0:raise ValueError('positive threshold')


def simulate_lifetime(s:LifetimeScenario):
    imbalance=0.;activity=1.;inventory=1.;trace=[];remix_times=[]
    for cy in range(1,s.cycles+1):
        activity*=1-s.felt_degradation_per_cycle
        inventory*=1-s.irreversible_inventory_loss_per_cycle
        imbalance+=(1-imbalance)*s.imbalance_growth_per_cycle
        pre_capacity=activity*inventory*(1-imbalance)
        mixed=False
        # When irreversible fade makes the target unreachable, repeated mixing
        # every cycle is not realistic. Only mix if some reversible imbalance
        # could physically account for the shortfall, with a minimum target.
        if s.remix_enabled and pre_capacity<s.remix_threshold_fraction and imbalance>.05:
            imbalance*=1-s.recovered_imbalance_fraction
            mixed=True;remix_times.append(cy)
        capacity=activity*inventory*(1-imbalance)
        trace.append({'cycle':cy,'pre_remix_capacity_fraction':pre_capacity,
                      'capacity_fraction':capacity,'felt_activity_fraction':activity,
                      'active_inventory_fraction':inventory,
                      'imbalance_fraction':imbalance,'remixed':mixed})
    return {'trace':trace,'remix_cycles':remix_times,
            'final_capacity_fraction':trace[-1]['capacity_fraction'],
            'final_felt_activity_fraction':activity,
            'final_active_inventory_fraction':inventory,
            'caveat':'illustrative phenomenological fade and idealized remix; not calibrated ageing'}
