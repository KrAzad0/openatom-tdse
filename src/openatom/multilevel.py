"""Six bound hydrogen states coupled to a prescribed z-directed drive (atomic units)."""
from math import factorial
import numpy as np
from scipy.integrate import quad
from scipy.special import eval_genlaguerre
from scipy.sparse import eye
from scipy.sparse.linalg import splu
from .torus import REDUCED_MASS, PROTON_MASS, HYDROGEN_MASS

LABELS=("1s","2s","2p0","3s","3p0","3d0")
NL=((1,0),(2,0),(2,1),(3,0),(3,1),(3,2))

def radial(n,l,r,mu=REDUCED_MASS):
    """Normalized real Coulomb radial function, with standard Laguerre phase."""
    if not isinstance(n,(int,np.integer)) or not isinstance(l,(int,np.integer)) or not 0<=l<n or not np.isfinite(mu) or mu<=0:
        raise ValueError("invalid hydrogen quantum numbers or reduced mass")
    rho=2*mu*np.asarray(r)/n
    return (2*mu/n)**1.5*np.sqrt(factorial(n-l-1)/(2*n*factorial(n+l)))*np.exp(-rho/2)*rho**l*eval_genlaguerre(n-l-1,2*l+1,rho)

def hydrogen_basis(mu=REDUCED_MASS):
    """Return excitation energies, z and radius matrices for LABELS."""
    if not np.isfinite(mu) or mu<=0:
        raise ValueError("invalid reduced mass")
    n=np.array([i for i,j in NL])
    e=-mu/(2*n*n)
    d=np.zeros((6,6))
    radius=d.copy()
    for i,(ni,li) in enumerate(NL):
        for j,(nj,lj) in enumerate(NL[:i+1]):
            if li==lj or abs(li-lj)==1:
                val=quad(lambda r: radial(ni,li,r,mu)*radial(nj,lj,r,mu)*r**3,0,np.inf,epsabs=1e-10)[0]
                if li==lj:
                    radius[i,j]=radius[j,i]=val
                else:
                    k=max(li,lj)
                    d[i,j]=d[j,i]=val*k/np.sqrt(4*k*k-1)
    return e-e[0],d,radius

class MultiStateTorus:
    """Nested symmetric splitting matching the six-state MATLAB method.

    H_int=diag(energies)-E*D; E is a signed drive, not a Maxwell solution.
    State rows index COM nodes; columns index internal states.
    """
    def __init__(self,grid,dt,energies=None,dipole=None):
        if not np.isfinite(dt) or dt<=0:
            raise ValueError("dt must be positive")
        if (energies is None)!=(dipole is None):
            raise ValueError("provide both energies and dipole or neither")
        if energies is None:
            energies,dipole,_=hydrogen_basis()
        e=np.asarray(energies)
        d=np.asarray(dipole)
        if e.ndim!=1 or len(e)<1 or d.shape!=(len(e),len(e)) or np.iscomplexobj(e) or np.iscomplexobj(d) or not np.isfinite(e).all() or not np.isfinite(d).all() or not np.allclose(d,d.T,rtol=0,atol=1e-13):
            raise ValueError("real energies and symmetric real dipole required")
        self.grid,self.dt,self.energies=grid,dt,e.copy()
        self.dipole=(d+d.T)/2
        self.phase=np.exp(-.5j*dt*e)
        self.lam,self.u=np.linalg.eigh(self.dipole)
        ident=eye(grid.kinetic.shape[0],dtype=complex,format="csc")
        self.right=(ident-1j*dt*grid.kinetic/4).tocsr()
        self.left=splu((ident+1j*dt*grid.kinetic/4).tocsc())

    def initial(self,sigma,momentum=0):
        psi=np.zeros((len(self.grid.xyz),len(self.energies)),complex)
        psi[:,0]=self.grid.packet(sigma,momentum)[:,0]
        return psi

    def step(self,state,field_midpoint):
        psi=np.asarray(state,complex)
        f=np.asarray(field_midpoint)
        n=len(self.grid.xyz)
        if psi.shape!=(n,len(self.energies)) or not np.isfinite(psi).all() or np.iscomplexobj(f) or not np.isfinite(f).all():
            raise ValueError("invalid state or field")
        f=np.broadcast_to(f,(n,))
        v=self.left.solve(self.right@psi)*self.phase
        v=((v@self.u)*np.exp(1j*self.dt*f[:,None]*self.lam))@self.u.T
        return self.left.solve(self.right@(v*self.phase))

    def observe(self,state,radius=None):
        psi=np.asarray(state,complex)
        if psi.shape!=(len(self.grid.xyz),len(self.energies)) or not np.isfinite(psi).all():
            raise ValueError("invalid state")
        dv=self.grid.dv
        prob=abs(psi)**2
        populations=dv*prob.sum(axis=0)
        norm=populations.sum()
        if norm<=0:
            raise ValueError("zero state")
        w=dv*prob.sum(axis=1)/norm
        com=w@self.grid.xyz
        rho=psi.conj().T@psi*dv/norm
        z=float(np.trace(rho@self.dipole).real)
        result=dict(populations=populations,norm=float(norm),com=com,relative_z=z,
                    electron_z=float(com[2]+PROTON_MASS/HYDROGEN_MASS*z),
                    proton_z=float(com[2]-z/HYDROGEN_MASS),
                    internal_energy=float(populations@self.energies/norm),
                    kinetic_energy=float(np.vdot(psi,self.grid.kinetic@psi).real*dv/norm))
        if radius is not None:
            result["mean_radius"]=float(np.trace(rho@radius).real)
        return result
