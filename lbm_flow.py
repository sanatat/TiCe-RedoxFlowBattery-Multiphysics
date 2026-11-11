"""Educational D2Q9 Lattice-Boltzmann BGK flow solver.
Periodic streamwise direction; no-slip horizontal walls and optional baffle obstacles.
Dimensionless *lattice units* only. Not a calibrated physical electrochemical CFD model.
"""
import numpy as np

C = np.array([[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]],dtype=int)
W = np.array([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36],float)
OPP = np.array([0,3,4,1,2,7,8,5,6])

def simulate_lbm(nx=100,ny=34,steps=2400,force=2e-6,tau=.85,with_baffles=False):
    solid=np.zeros((ny,nx),bool);solid[0]=True;solid[-1]=True
    if with_baffles:
        solid[1:ny//3,nx//3:nx//3+5]=True
        solid[2*ny//3:ny-1,2*nx//3:2*nx//3+5]=True
    f=np.broadcast_to(W[:,None,None],(9,ny,nx)).copy()
    Fx=float(force)
    for step in range(steps):
        rho=np.sum(f,axis=0)
        momx=np.tensordot(C[:,0],f,axes=(0,0))
        momy=np.tensordot(C[:,1],f,axes=(0,0))
        ux=(momx+Fx*.5)/rho;uy=momy/rho
        ux[solid]=0;uy[solid]=0
        u2=ux*ux+uy*uy
        dot=C[:,0,None,None]*ux+C[:,1,None,None]*uy
        feq=W[:,None,None]*rho[None,:,:]*(1+3*dot+4.5*dot**2-1.5*u2)
        # Guo body-force scheme for horizontal pressure-gradient equivalence
        guo=W[:,None,None]*(1-.5/tau)*(3*(C[:,0,None,None]-ux)+9*dot*C[:,0,None,None])*Fx
        post=f+(feq-f)/tau+guo
        # on-site no-slip bounce-back at obstacles; streamed in subsequent step
        post[:,solid]=f[OPP][:,solid]
        f=np.stack([np.roll(post[i],shift=(cy,cx),axis=(0,1)) for i,(cx,cy) in enumerate(C)])
    rho=np.sum(f,axis=0)
    ux=(np.tensordot(C[:,0],f,axes=(0,0))+.5*Fx)/rho
    uy=np.tensordot(C[:,1],f,axes=(0,0))/rho
    ux[solid]=0.;uy[solid]=0.
    return {'ux':ux,'uy':uy,'solid':solid,'mean_ux':float(ux[~solid].mean()),
            'density_range':[float(rho.min()),float(rho.max())],
            'tau':tau,'nu_lattice':(tau-.5)/3,'force_lattice':force,
            'steps':steps,'with_baffles':with_baffles}

if __name__=='__main__':
    for k in (False, True):
        import time
        st=time.time();r=simulate_lbm(steps=2000,with_baffles=k)
        print(k,r['mean_ux'],r['density_range'],'elapsed',time.time()-st)
