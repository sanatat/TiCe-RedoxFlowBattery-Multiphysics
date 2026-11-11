"""Virtual semi-infinite CV via implicit finite-difference diffusion + BV Robin.
Uncompensated resistance and double-layer capacitance omitted; only a
conceptual virtual 3-electrode experiment, not a fit to article CV curves.
"""
import numpy as np
from scipy.linalg import solve_banded
F=96485.33212;R=8.314462618

def simulate(D=1.78e-8,k0=8.1e-6,conc=100.,scan=.005,T=298.15,
             xmax=.02,ny=240,dt=.25):
    dx=xmax/ny
    t_leg=.5/scan
    ts=np.arange(0,2*t_leg+dt/2,dt)
    e=np.where(ts<=t_leg,-.25+scan*ts,.25-scan*(ts-t_leg))
    c=np.zeros(ny) # oxidized concentration, reduced bulk initial
    ratio=D*dt/dx**2
    ab=np.zeros((3,ny));ab[0,1:]=-ratio;ab[2,:-1]=-ratio
    records=[]
    for ix,potential in enumerate(e):
        z=np.clip(.5*F*potential/(R*T),-25,25)
        a=np.exp(z);b=np.exp(-z)
        ksum=k0*(a+b)
        # Robin flux J = k0 (a*(C-cs) - b*cs); half-cell centre correction
        eff=1/(1+ksum*dx/(2*D))
        source=k0*conc*a*eff
        sink=ksum*eff
        diag=np.full(ny,1+2*ratio)
        diag[0]=1+ratio+dt*sink/dx
        diag[-1]=1+ratio # remote zero flux; front must not reach there
        ab[1]=diag
        rhs=c.copy();rhs[0]+=dt*source/dx
        c=solve_banded((1,1),ab,rhs)
        if ix%4==0:
            flux=source-sink*c[0]
            records.append((float(potential),float(F*flux)))
    return np.array(records)
