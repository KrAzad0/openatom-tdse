import numpy as np
from openatom import CrankNicolson, gaussian_pulse

def test_norm_time_dependent():
    x = np.linspace(-10, 10, 401)[1:-1]
    p = CrankNicolson(x, 0.01, lambda x,t: x*x/2 + x*gaussian_pulse(t, .1, 2, .5, .2))
    psi = np.exp(-x*x/2).astype(complex)
    psi /= np.sqrt(p.norm(psi))
    for k in range(100):
        psi = p.step(psi, k*p.dt)
    np.testing.assert_allclose(p.norm(psi), 1, atol=1e-12)

def test_discrete_free_eigenmode():
    n, length, dt = 199, 4., .02
    dx = length/(n+1)
    x = np.arange(1,n+1)*dx
    p = CrankNicolson(x, dt, lambda x,t: 0.)
    psi = np.sin(np.pi*x/length)
    energy = (1-np.cos(np.pi/(n+1)))/dx**2
    factor = (1-.5j*dt*energy)/(1+.5j*dt*energy)
    np.testing.assert_allclose(p.step(psi, 0), factor*psi, atol=1e-12)
