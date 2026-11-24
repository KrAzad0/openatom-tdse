"""Illustrative 1D soft-core atom; not a quantitative 3D hydrogen model."""
import numpy as np
from scipy.linalg import eigh_tridiagonal
from openatom import CrankNicolson, gaussian_pulse

x = np.linspace(-40, 40, 1601)[1:-1]
dx = x[1]-x[0]
v = -1/np.sqrt(x*x+1)
e, u = eigh_tridiagonal(1/dx**2+v, np.full(x.size-1,-.5/dx**2), select="i", select_range=(0,0))
psi = u[:,0].astype(complex)/np.sqrt(dx)
p = CrankNicolson(x, .02, lambda x,t: v+x*gaussian_pulse(t,.05,2.,2.,.5))
for k in range(200):
    psi = p.step(psi, k*p.dt)
print("Initial energy (Ha):", e[0])
print("Final norm:", p.norm(psi))
print("Final dipole <x> (bohr):", dx*np.sum(x*abs(psi)**2))
