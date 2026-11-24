"""Radial Coulomb eigenstates, in atomic units."""
import numpy as np
from scipy.linalg import eigh_tridiagonal

def hydrogen_radial(l=0, Z=1.0, r_max=80.0, points=2000, states=3):
    """Return interior radii, energies, and normalized reduced radial states.

    Dirichlet boundaries at r=0 and r=r_max. Columns of u satisfy
    integral |u|^2 dr ~= 1; principal quantum number n=l+1+j
    for the bound eigenstates. Positive box states are not bound states.
    """
    if not isinstance(l, (int, np.integer)) or l < 0:
        raise ValueError("l must be a nonnegative integer")
    if not isinstance(points, (int, np.integer)) or not isinstance(states, (int, np.integer)):
        raise ValueError("points and states must be integers")
    if not np.isfinite([Z, r_max]).all() or Z <= 0 or r_max <= 0 or points < 2 or not 1 <= states <= points:
        raise ValueError("invalid grid, nuclear charge or state count")
    h = r_max / (points + 1)
    r = h * np.arange(1, points + 1)
    diagonal = 1 / h**2 + l * (l + 1) / (2 * r**2) - Z / r
    off = np.full(points - 1, -0.5 / h**2)
    energies, u = eigh_tridiagonal(diagonal, off, select="i", select_range=(0, states - 1))
    return r, energies, u / np.sqrt(h)
