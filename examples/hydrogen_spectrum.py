"""Run after pip install -e ."""
from openatom import hydrogen_radial

r, energies, u = hydrogen_radial()
print("n  numerical (Ha)  analytical (Ha)  absolute error")
for n, e in enumerate(energies, 1):
    exact = -0.5/n**2
    print(f"{n}  {e:.8f}  {exact:.8f}  {abs(e-exact):.3e}")
