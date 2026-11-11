"""Conservative finite-volume 2-D advection/diffusion/first-order wall reaction.

Research training simulator. Simplified flow channel, NOT the paper's unresolved
porous carbon-felt flow field. Units: SI; velocity either analytical plane
Poiseuille or lattice-Boltzmann field mapped onto an *assumed* channel.
"""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

@dataclass(frozen=True)
class Channel:
    length_m: float = 0.03       # assumed physical length, not documented channel
    height_m: float = 0.001     # assumed fluid-channel gap, not felt thickness
    width_m: float = 0.03       # equal to geometric cell width, assumed flow path
    flow_m3_s: float = 20e-6/60 # paper flow in each half cell
    mu_Pa_s: float = 2.76e-3   # SI S7 0.5M, room temperature
    rho_kg_m3: float = 1100.0  # illustrative density, not reported
    efficiency: float = 0.55  # illustrative pump efficiency

    @property
    def u_mean(self):return self.flow_m3_s/(self.width_m*self.height_m)
    @property
    def reynolds(self):return self.rho_kg_m3*self.u_mean*(2*self.height_m)/self.mu_Pa_s
    @property
    def pressure_drop_Pa(self):return 12*self.mu_Pa_s*self.length_m*self.u_mean/self.height_m**2
    @property
    def pump_power_two_cells_W(self):return 2*self.pressure_drop_Pa*self.flow_m3_s/self.efficiency


def poiseuille_velocity(ch:Channel,ny:int,nx:int):
    """Cell-centred incompressible analytical velocity, normalized to exact Q."""
    y=(np.arange(ny)+.5)/ny
    profile=6*y*(1-y)
    profile=profile/profile.mean()
    ux=np.broadcast_to((ch.u_mean*profile)[:,None],(ny,nx)).copy()
    uy=np.zeros_like(ux)
    solid=np.zeros_like(ux,dtype=bool)
    return ux,uy,solid


def lbm_velocity(ch:Channel,nx=84,ny=28,baffles=False,steps=1600):
    """A D2Q9 field scaled to the assumed physical channel's prescribed inlet flow.

    LBM runs on a periodic segment in lattice units. Mapping x/y to the
    user-specified length/height is an illustrative coordinate transform,
    *not* a calibrated hydrodynamics solution when dx_physical != dy_physical.
    """
    from lbm_flow import simulate_lbm
    out=simulate_lbm(nx=nx,ny=ny,steps=steps,with_baffles=baffles)
    ux=out['ux'].copy();uy=out['uy'].copy();solid=out['solid']
    # Scale to match actual incoming volumetric flow across fluid inlet.
    inlet_mean=ux[~solid[:,0],0].sum()/ny
    if inlet_mean<=0:raise RuntimeError('LBM inlet velocity not positive')
    scale=ch.u_mean/inlet_mean
    ux*=scale;uy*=scale*((ch.height_m/ny)/(ch.length_m/nx)) # preserve lattice aspect ratio in anisotropic physical mapping
    ux[solid]=0;uy[solid]=0
    return ux,uy,solid,out


def channel_transport(ch:Channel,D_m2_s:float,cin_mol_m3:float,
                      k_wall_m_s:float=1.0,nx:int=96,ny:int=32,
                      flow_field=None):
    """Steady FVM with inlet concentration, zero diffusive outlet flux,
    no-through-flux walls/obstacles and Robin consumption at bottom electrode.

    Returns a strict global scalar balance by construction. A nonzero reported
    balance from the PDE solve indicates a numerical error. Outflow and inlet
    are evaluated with the same discrete fluxes as the linear system.

    Reaction flux j(x)=k_wall*c_surface, effective k=(1/k+dy/2D)^-1.
    Electroactive bottom wall is only a *surrogate* for porous felt.
    """
    if min(D_m2_s,cin_mol_m3,k_wall_m_s)<=0:raise ValueError('positive D,c,k required')
    if flow_field is None: ux,uy,solid=poiseuille_velocity(ch,ny,nx)
    else:
        ux,uy,solid=flow_field
        ny,nx=ux.shape
    dx=ch.length_m/nx;dy=ch.height_m/ny
    idx=np.full((ny,nx),-1,int);idx[~solid]=np.arange(np.sum(~solid));n=idx.max()+1
    ar=[];ac=[];av=[];rhs=np.zeros(n)
    def add(i,j,val):ar.append(i);ac.append(j);av.append(val)
    def coupled(i,j,uarea,cond):
        """Flux i->j = adv_upwind + conductance*(ci-cj)."""
        if uarea>=0:
            add(i,i,uarea+cond);add(i,j,-cond)
            add(j,i,-uarea-cond);add(j,j,cond)
        else:
            add(i,i,cond);add(i,j,uarea-cond)
            add(j,i,-cond);add(j,j,-uarea+cond)
    # Iterate internal vertical and horizontal faces, no flux through solids.
    for y in range(ny):
        for x in range(nx):
            i=idx[y,x]
            if i<0:continue
            if x+1<nx and idx[y,x+1]>=0:
                j=idx[y,x+1]
                uarea=.5*(ux[y,x]+ux[y,x+1])*dy
                coupled(i,j,uarea,D_m2_s*dy/dx)
            if y+1<ny and idx[y+1,x]>=0:
                j=idx[y+1,x]
                varea=.5*(uy[y,x]+uy[y+1,x])*dx
                coupled(i,j,varea,D_m2_s*dx/dy)
    # inlet/outlet boundaries
    q_in_coeff=0;diff_in_coeff=0;q_out_coeffs=[]
    for y in range(ny):
        i=idx[y,0]
        if i>=0:
            adv=ux[y,0]*dy
            cond=2*D_m2_s*dy/dx
            if adv>=0:
                add(i,i,cond);rhs[i]+=(adv+cond)*cin_mol_m3
                q_in_coeff+=adv;diff_in_coeff+=cond
            else:
                add(i,i,cond-adv);rhs[i]+=cond*cin_mol_m3
                diff_in_coeff+=cond
        j=idx[y,-1]
        if j>=0:
            adv=ux[y,-1]*dy
            if adv>=0:
                add(j,j,adv);q_out_coeffs.append((j,adv))
            else:
                rhs[j]+=-adv*cin_mol_m3 # prescribed backflow concentration
    # reaction along bottom (only fluid cells), flux wall = k_eff * c0 * dx
    k_eff=1/(1/k_wall_m_s+dy/(2*D_m2_s))
    wall_coeffs=[]
    if not solid[0].all():
        for x in range(nx):
            i=idx[0,x]
            if i>=0:
                add(i,i,k_eff*dx);wall_coeffs.append((i,k_eff*dx))
    else:
        # Bottom LBM bounce-back solid row, reactive wall adjacent to solid.
        for x in range(nx):
            for y in range(ny):
                if not solid[y,x] and (y==0 or solid[y-1,x]):
                    # only bottom wall, not side of internal baffles
                    if y<=1:
                        i=idx[y,x];add(i,i,k_eff*dx);wall_coeffs.append((i,k_eff*dx))
                    break
    A=coo_matrix((av,(ar,ac)),shape=(n,n)).tocsr()
    sol=spsolve(A,rhs)
    if not np.all(np.isfinite(sol)) or sol.min() < -1e-7*cin_mol_m3:
        raise ArithmeticError('transport nonphysical concentration')
    c=np.full((ny,nx),np.nan);c[~solid]=sol
    # same boundary flux definitions as used above
    inlet=0.0;outlet=0.0
    for y in range(ny):
        i=idx[y,0]
        if i>=0:
            adv=ux[y,0]*dy;cond=2*D_m2_s*dy/dx
            inlet+=(adv*cin_mol_m3 if adv>=0 else adv*sol[i])+cond*(cin_mol_m3-sol[i])
        j=idx[y,-1]
        if j>=0:
            adv=ux[y,-1]*dy
            outlet+=adv*(sol[j] if adv>=0 else cin_mol_m3)
    uptake=sum(s*sol[i] for i,s in wall_coeffs)
    error=(inlet-outlet-uptake)/max(abs(inlet),1e-30)
    return {'c':c,'ux':ux,'uy':uy,'solid':solid,'x_m':(np.arange(nx)+.5)*dx,
            'y_m':(np.arange(ny)+.5)*dy,'in_mol_s':inlet*ch.width_m,
            'out_mol_s':outlet*ch.width_m,'reacted_mol_s':uptake*ch.width_m,
            'uptake_kmol_s':uptake*ch.width_m*1e-3,
            'mass_balance_relative_error':float(error),
            'k_apparent_m_s':uptake/(ch.length_m*cin_mol_m3),
            'j_equiv_mA_cm2':96485.33212*uptake*ch.width_m/(0.0009*10),
            'minimum_c_mol_m3':float(np.nanmin(c)),
            'delta_p_Pa':ch.pressure_drop_Pa,'pump_power_two_cells_W':ch.pump_power_two_cells_W,
            'reynolds':ch.reynolds,'nx':nx,'ny':ny}
