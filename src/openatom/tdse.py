"""One-dimensional midpoint Crank-Nicolson propagation in atomic units."""
import numpy as np
from scipy.linalg import solve_banded

def gaussian_pulse(t, amplitude, omega, center, sigma, phase=0.0):
    """Electric field with Gaussian envelope; sigma is amplitude standard deviation."""
    if not np.isfinite([amplitude, omega, center, sigma, phase]).all() or sigma <= 0:
        raise ValueError("invalid pulse parameters")
    tau = np.asarray(t) - center
    return amplitude * np.exp(-0.5 * (tau / sigma)**2) * np.cos(omega * tau + phase)

class CrankNicolson:
    """Propagate on an equally spaced interior grid with zero exterior endpoints.

    potential(x,t) returns a real scalar or a real array of shape x.shape.
    No absorber is included; reflections at finite-box boundaries are possible.
    """
    def __init__(self, x, dt, potential):
        self.x = np.asarray(x, dtype=float).copy()
        self.dt = float(dt)
        self.potential = potential
        if self.x.ndim != 1 or self.x.size < 2 or not np.isfinite(self.x).all():
            raise ValueError("invalid spatial grid")
        dx = np.diff(self.x)
        if dx[0] <= 0 or not np.allclose(dx, dx[0], rtol=1e-9, atol=0):
            raise ValueError("grid must be uniform and increasing")
        if not np.isfinite(self.dt) or self.dt <= 0:
            raise ValueError("dt must be positive")
        self.dx = dx[0]
        self.off = -0.5 / self.dx**2

    def norm(self, psi):
        return float(self.dx * np.sum(np.abs(psi)**2))

    def step(self, psi, t):
        psi = np.asarray(psi, dtype=complex)
        if psi.shape != self.x.shape or not np.isfinite(psi).all() or not np.isfinite(t):
            raise ValueError("invalid state or time")
        v = np.asarray(self.potential(self.x, t + self.dt / 2))
        if np.iscomplexobj(v) or not np.isfinite(v).all():
            raise ValueError("potential must be finite and real")
        v = np.broadcast_to(v, self.x.shape)
        d = 1 / self.dx**2 + v
        a = 0.5j * self.dt
        rhs = (1 - a*d) * psi
        rhs[1:] -= a*self.off*psi[:-1]
        rhs[:-1] -= a*self.off*psi[1:]
        bands = np.zeros((3, self.x.size), dtype=complex)
        bands[0, 1:] = a*self.off
        bands[1] = 1 + a*d
        bands[2, :-1] = a*self.off
        return solve_banded((1, 1), bands, rhs)
