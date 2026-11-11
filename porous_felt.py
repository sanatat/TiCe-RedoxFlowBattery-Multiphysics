"""Analytical and FVM check of 1-D porous carbon-felt diffusion/reaction.

Assumed first-order irreversible, effective D and electroactive area density;
not the article's validated carbon felt microstructure or reversible BV model.
This submodel investigates reactive penetration/utilization, not full battery kinetics.
"""
from dataclasses import dataclass
import numpy as np
from electrochem_cell import published_apparent_D,published_apparent_k

@dataclass(frozen=True)
class Felt:
    thickness_m:float=.005  # SI: measured full-cell felt thickness
    porosity:float=.85       # illustrative, not in SI
    tortuosity:float=2.5     # illustrative, not in SI
    active_area_m2_m3:float=25000. # illustrative, not in SI


def effectiveness(species='Ti3',f=Felt(),D_scale=1.0,c_surface=250.):
    D_eff=published_apparent_D(species)*D_scale*f.porosity/f.tortuosity
    k_surface=published_apparent_k(species)
    rate_vol=k_surface*f.active_area_m2_m3
    phi=f.thickness_m*np.sqrt(rate_vol/D_eff)
    eta=np.tanh(phi)/phi if phi>1e-12 else 1.
    k_eff=np.sqrt(D_eff*rate_vol)*np.tanh(phi)
    z=np.linspace(0,f.thickness_m,140)
    # Overflow-safe cosh ratio for phi possibly very high.
    zz=np.exp(-phi*z/f.thickness_m)
    profile=c_surface*(zz+np.exp(-2*phi)*np.exp(phi*z/f.thickness_m))/(1+np.exp(-2*phi))
    return {'D_eff_m2_s':D_eff,'k_s_m_s':k_surface,'k_vol_s_inv':rate_vol,
            'Thiele':phi,'effectiveness':eta,'k_wall_m_s':k_eff,
            'penetration_depth_m':np.sqrt(D_eff/rate_vol),
            'z_m':z,'concentration_mol_m3':profile,
            'flux_mol_m2_s':k_eff*c_surface}


def finite_difference_flux(species='Ti3',f=Felt(),D_scale=1.0,c_surface=250.,n=180):
    """Independent finite-volume slab solution: Dirichlet at wetted inlet,
    no-flux collector at back of felt. Compare to analytical solution.
    """
    from scipy.sparse import diags
    from scipy.sparse.linalg import spsolve
    d=effectiveness(species,f,D_scale,c_surface)
    D=d['D_eff_m2_s'];rate=d['k_vol_s_inv'];dx=f.thickness_m/n
    coeff=D/dx**2
    diagonals=np.full(n,2*coeff+rate)
    diagonals[0]=3*coeff+rate
    diagonals[-1]=coeff+rate
    off=np.full(n-1,-coeff)
    A=diags([off,diagonals,off],[-1,0,1],format='csr')
    rhs=np.zeros(n);rhs[0]=2*coeff*c_surface
    c=spsolve(A,rhs)
    flux=2*D*(c_surface-c[0])/dx
    return {'flux_mol_m2_s':flux,'analytical_flux':d['flux_mol_m2_s'],
            'relative_error':abs(flux-d['flux_mol_m2_s'])/d['flux_mol_m2_s'],
            'minimum_concentration':c.min()}
