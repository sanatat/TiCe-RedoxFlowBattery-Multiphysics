"""Conservative two-tank membrane cross-over and hypothetical water migration.

Ce effective permeability from SI S4. Ion charge is effective complex charge,
not bare Ce(III)/Ti(IV) formal valence: actual distributions are unknown.
Electric migration uses a *linearized dilute Nernst-Planck scenario* and a
prescribed membrane potential drop. Electroneutrality, supporting-ion motion,
electro-osmotic drag and membrane partitioning are NOT solved. Interpret as
sensitivity to missing mechanisms rather than a prediction of actual flux.
"""
from dataclasses import dataclass
import numpy as np
from electrochem_cell import Battery,F,R

@dataclass(frozen=True)
class MembraneScenario:
    volume_each_m3:float=20e-6
    area_m2:float=9e-4
    Ce_permeability_m_s:float=4.54e-9 # SI S4 measured in H-cell
    Ti_permeability_m_s:float=2e-10 # scenario; below detection != zero
    z_Ce_effective:float=0.0 # hypothesis: complex charge uncertain
    z_Ti_effective:float=0.0
    potential_drop_V:float=0.0 # prescribed, not measured voltage across N115
    water_transfer_m3_s:float=0.0 # positive moves tank A -> tank B
    temperature_K:float=298.15

    def __post_init__(self):
        if min(self.volume_each_m3,self.area_m2,self.temperature_K)<=0:raise ValueError('positive geometry/T')
        if min(self.Ce_permeability_m_s,self.Ti_permeability_m_s)<0:raise ValueError('nonnegative P')


def simulate_membrane(s:MembraneScenario,hours=120.,dt_s=60.,
                      initial_concentrations_mol_m3=((500.,0.),(0.,500.))):
    """Track total Ce and Ti moles, not individual redox states.

    Fixed-gradient cross-membrane flux approximation with substepping for
    positive inventories; when water shifts, concentrations=moles/volumes.
    Chemical electron/charge conservation is outside this model's scope.
    """
    initial=np.asarray(initial_concentrations_mol_m3,dtype=float)
    if initial.shape!=(2,2) or np.min(initial)<0:raise ValueError('2 tanks x [Ce,Ti] nonnegative concentrations')
    vols=np.array([s.volume_each_m3,s.volume_each_m3],float)
    n=initial*vols[:,None]
    init_inventory=n.sum(axis=0).copy();init_volume=vols.sum()
    timeline=[]
    steps=max(1,int(np.ceil(hours*3600/dt_s)))
    dt=hours*3600/steps
    p=np.array([s.Ce_permeability_m_s,s.Ti_permeability_m_s])
    z=np.array([s.z_Ce_effective,s.z_Ti_effective])
    for k in range(steps+1):
        if k==0 or k==steps or k%max(1,steps//120)==0:
            c=n/vols[:,None]
            timeline.append({'time_h':k*dt/3600,'Ce_A_M':c[0,0]/1000,'Ce_B_M':c[1,0]/1000,
                             'Ti_A_M':c[0,1]/1000,'Ti_B_M':c[1,1]/1000,
                             'volume_A_mL':vols[0]*1e6,'volume_B_mL':vols[1]*1e6})
        if k==steps:break
        c=n/vols[:,None]
        flux=p*(c[0]-c[1]) + p*z*F*s.potential_drop_V/(R*s.temperature_K)*.5*(c[0]+c[1])
        # 2x2 inventories, equal/opposite fluxes conserve EACH element.
        dn=flux*s.area_m2*dt
        # Bounded to avoid creating negative inventories at large fields.
        dn=np.minimum(np.maximum(dn,-n[1]),n[0])
        n[0]-=dn;n[1]+=dn
        dv=s.water_transfer_m3_s*dt
        if min(vols[0]-dv,vols[1]+dv)<=0:raise ValueError('water scenario exhausts a reservoir')
        vols[0]-=dv;vols[1]+=dv
    return {'trace':timeline,'initial_Ce_mol':init_inventory[0],
            'final_Ce_mol':n[:,0].sum(),'initial_Ti_mol':init_inventory[1],
            'final_Ti_mol':n[:,1].sum(),'initial_volume_m3':init_volume,
            'final_volume_m3':vols.sum(),'final_tank_moles':n,'final_volumes_m3':vols,
            'scenario':s.__dict__}
