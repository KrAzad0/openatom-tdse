"""Experimental two-state hydrogen COM dynamics in a solid torus (atomic units)."""
import numpy as np
from scipy.sparse import coo_matrix, eye
from scipy.sparse.linalg import splu

PROTON_MASS = 1836.15267343
HYDROGEN_MASS = 1.0 + PROTON_MASS
REDUCED_MASS = PROTON_MASS / HYDROGEN_MASS
ATOMIC_TIME_AS = 24.18884326505
GAP = 3 * REDUCED_MASS / 8
DIPOLE = 256 / (243 * np.sqrt(2) * REDUCED_MASS)

class TorusGrid:
    """Cartesian interior nodes with a zero-exterior seven-point stencil."""
    def __init__(self, major, minor, spacing, mass=HYDROGEN_MASS):
        if not np.isfinite([major, minor, spacing, mass]).all() or not major > minor > 0 or min(spacing,mass) <= 0:
            raise ValueError("require major > minor > 0, positive spacing and mass")
        self.major, self.minor, self.h = major, minor, spacing
        self.dv = spacing**3
        nxy = int(np.ceil((major+minor)/spacing))
        nz = int(np.ceil(minor/spacing))
        axes = [np.arange(-nxy,nxy+1), np.arange(-nxy,nxy+1), np.arange(-nz,nz+1)]
        X,Y,Z = np.meshgrid(*axes,indexing="ij")
        mask = (np.hypot(X*spacing,Y*spacing)-major)**2+(Z*spacing)**2 < minor**2
        sites = np.argwhere(mask)
        self.xyz = np.column_stack((X[mask],Y[mask],Z[mask]))*spacing
        n = len(sites)
        if n < 2:
            raise ValueError("grid too coarse")
        lookup = np.full(mask.shape,-1,dtype=int)
        lookup[mask] = np.arange(n)
        rows, cols = [np.arange(n)], [np.arange(n)]
        vals = [np.full(n,3/(mass*spacing**2))]
        for axis in range(3):
            for sign in (-1,1):
                other = sites.copy()
                other[:,axis] += sign
                valid = (other[:,axis]>=0)&(other[:,axis]<mask.shape[axis])
                src = np.flatnonzero(valid)
                dst = lookup[tuple(other[valid].T)]
                keep = dst>=0
                rows.append(src[keep]); cols.append(dst[keep])
                vals.append(np.full(keep.sum(),-0.5/(mass*spacing**2)))
        self.kinetic = coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(n,n)).tocsr()

    @property
    def volume(self):
        return len(self.xyz)*self.dv

    @property
    def exact_volume(self):
        return 2*np.pi**2*self.major*self.minor**2

    def profile(self, angular=0.0, vertical=0.0):
        if not np.isfinite([angular,vertical]).all() or abs(angular)+abs(vertical)>=1:
            raise ValueError("require |angular| + |vertical| < 1")
        x,y,z = self.xyz.T
        return 1+angular*x/np.hypot(x,y)+vertical*z/self.minor

    def packet(self, sigma, momentum=0.0):
        """Density standard deviations and initial tangential phase momentum."""
        sigma=np.asarray(sigma,float)
        if sigma.shape!=(3,) or not np.isfinite(sigma).all() or np.any(sigma<=0) or not np.isfinite(momentum):
            raise ValueError("invalid packet parameters")
        displacement=self.xyz-np.array([self.major,0,0])
        exponent=-0.25*np.sum((displacement/sigma)**2,axis=1)
        psi=np.exp(exponent-exponent.max()+1j*momentum*self.xyz[:,1])
        state=np.zeros((len(psi),2),complex)
        state[:,0]=psi/np.sqrt(self.dv*np.sum(abs(psi)**2))
        return state

class TwoStateTorus:
    """Strang split: kinetic half Cayley steps, exact local midpoint two-state step.

    Convention: H_int = [[0,-d E],[-d E,gap]], as in the original
    two-state model. E is the signed drive parameter in this convention.
    """
    def __init__(self, grid, dt, gap=GAP, dipole=DIPOLE):
        if not np.isfinite([dt,gap,dipole]).all() or dt<=0 or gap<=0:
            raise ValueError("invalid timestep or internal parameters")
        self.grid,self.dt,self.gap,self.dipole=grid,dt,gap,dipole
        ident=eye(grid.kinetic.shape[0],format="csc",dtype=complex)
        self._right=(ident-1j*dt*grid.kinetic/4).tocsr()
        self._left=splu((ident+1j*dt*grid.kinetic/4).tocsc())

    def populations(self,state):
        return self.grid.dv*np.sum(abs(np.asarray(state))**2,axis=0)

    def step(self,state,field_midpoint):
        """Advance dt; field_midpoint is a real scalar or a field at each node."""
        state=np.asarray(state,complex)
        n=len(self.grid.xyz)
        field=np.asarray(field_midpoint)
        if state.shape!=(n,2) or not np.isfinite(state).all() or np.iscomplexobj(field) or not np.isfinite(field).all():
            raise ValueError("invalid state or electric field")
        g=-self.dipole*np.broadcast_to(field,(n,))
        a=self.gap/2
        omega=np.hypot(a,g)
        c=np.cos(self.dt*omega)
        s=-1j*self.dt*np.sinc(self.dt*omega/np.pi)
        v=self._left.solve(self._right@state)
        w=np.empty_like(v)
        w[:,0]=c*v[:,0]+s*(-a*v[:,0]+g*v[:,1])
        w[:,1]=c*v[:,1]+s*(g*v[:,0]+a*v[:,1])
        w*=np.exp(-1j*self.dt*a)
        return self._left.solve(self._right@w)
