"""Conditional mixed-acid *scenarios*, not verified Ti/Ce speciation.

SI S7/S8 quantify 2 M MSA-based electrolyte only.  Parameters describing
sulfate addition, metal complexation, and solubility are USER-ASSUMED.  They
must be replaced by experimental equilibrium constants / measurements before
any chemical precipitation or acid-optimum claim can be made.
"""
from dataclasses import dataclass
import math
from electrochem_cell import electrolyte_props


@dataclass(frozen=True)
class MixedAcidHypothesis:
    sulfate_fraction: float = 0.0  # fraction of nominal 2 M acid replaced by H2SO4
    sulfate_to_msa_conductivity_factor: float = 1.15  # hypothesis, NOT literature measurement
    sulfate_to_msa_viscosity_factor: float = 0.90  # hypothesis, NOT literature measurement
    conditional_beta_Ce_M_inv: float = 1.0  # hypothetical Ce/sulfate association
    conditional_beta_Ti_M_inv: float = 3.0  # hypothetical Ti/sulfate association
    intrinsic_Ce_solubility_M: float = 0.30  # arbitrary illustrative threshold
    intrinsic_Ti_solubility_M: float = 0.35  # arbitrary illustrative threshold

    def __post_init__(self):
        if not 0 <= self.sulfate_fraction <= 1:
            raise ValueError('sulfate_fraction must lie in [0,1]')
        if min(self.sulfate_to_msa_conductivity_factor,
               self.sulfate_to_msa_viscosity_factor,
               self.intrinsic_Ce_solubility_M,
               self.intrinsic_Ti_solubility_M) <= 0:
            raise ValueError('positive physical scenario parameters required')
        if min(self.conditional_beta_Ce_M_inv,self.conditional_beta_Ti_M_inv)<0:
            raise ValueError('association constants cannot be negative')


def conditional_partition(total_M, ligand_M, beta_M_inv):
    """Simple 1:1 conditional binding, [ML]=beta*[Mfree]*[Lfree].

    L_free is PROVIDED (not a solved sulfate acid-base equilibrium).
    This is a transparent sensitivity model, not actual Ti/Ce-MSA speciation.
    """
    if min(total_M,ligand_M,beta_M_inv)<0:raise ValueError('nonnegative inputs')
    free=total_M/(1+beta_M_inv*ligand_M)
    return {'free_M':free,'bound_M':total_M-free}


def mixed_acid_scenario(h:MixedAcidHypothesis,total_each_M=0.5,T_C=25.0):
    """Anchor at measured 2 M MSA properties, vary unknown sulfate effects."""
    if total_each_M <= 0:raise ValueError('positive metal loading')
    conductivity_mS_cm,mu_cP=electrolyte_props(total_each_M,T_C)
    conductivity_mS_cm *=10  # S/m -> mS/cm
    mu_cP*=1000 # Pa s -> cP
    f=h.sulfate_fraction
    sigma=conductivity_mS_cm*(1+(h.sulfate_to_msa_conductivity_factor-1)*f)
    mu=mu_cP*h.sulfate_to_msa_viscosity_factor**f
    sulfate_nominal_M=2.0*f # NOT free sulfate or proton activity
    ce=conditional_partition(total_each_M,sulfate_nominal_M,h.conditional_beta_Ce_M_inv)
    ti=conditional_partition(total_each_M,sulfate_nominal_M,h.conditional_beta_Ti_M_inv)
    # The following are *conditional screening ratios*, not saturation indices
    # calculated from Ti/Ce thermodynamic activity/solid-phase equilibrium.
    ce_capacity=h.intrinsic_Ce_solubility_M*(1+h.conditional_beta_Ce_M_inv*sulfate_nominal_M)
    ti_capacity=h.intrinsic_Ti_solubility_M*(1+h.conditional_beta_Ti_M_inv*sulfate_nominal_M)
    return {'sulfate_fraction':f,'MSA_nominal_M':2*(1-f),
            'H2SO4_nominal_M':sulfate_nominal_M,
            'conductivity_mS_cm_scenario':sigma,'viscosity_cP_scenario':mu,
            'Ce_free_M_conditional':ce['free_M'],'Ce_bound_M_conditional':ce['bound_M'],
            'Ti_free_M_conditional':ti['free_M'],'Ti_bound_M_conditional':ti['bound_M'],
            'Ce_loading_to_assumed_capacity_ratio':total_each_M/ce_capacity,
            'Ti_loading_to_assumed_capacity_ratio':total_each_M/ti_capacity,
            'thermodynamically_validated':False,
            'interpretation':'hypothetical screening only; NOT measured precipitation risk'}
