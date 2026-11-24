# Numerical theory

All equations use atomic units (electron mass, elementary charge, reduced Planck constant and Coulomb constant set to one). Length is in bohr, time in atomic time units, energy in hartree.

## Radial Coulomb problem

With u(r)=r R(r), the fixed-nucleus radial Hamiltonian is
H = -0.5 d²/dr² + l(l+1)/(2r²) - Z/r.
We impose u(0)=u(r_max)=0 and use centered second differences on the interior grid, giving a real symmetric tridiagonal matrix with off-diagonal -1/(2h²). The diagonal is 1/h² plus the effective potential. The analytic infinite-domain bound energies are -Z²/(2n²), n>=l+1. Normalize with h sum |u|²=1.

Check both grid-spacing and domain-size convergence. High positive-energy box states are not continuum-normalized states.

## Time propagation

For i dψ/dt = H(t)ψ, define H at the temporal midpoint and solve
(I+i dt H/2)ψ_next = (I-i dt H/2)ψ.
For real potentials the discrete Hamiltonian is Hermitian and this Cayley transform is unitary in exact arithmetic. The midpoint method is second-order in time for sufficiently smooth H. Verify convergence independently; norm conservation alone does not establish physical accuracy.

The soft-core example uses V(x)=-1/sqrt(x²+1) and length-gauge interaction +x E(t). Its Gaussian electric field is illustrative only, not a calibrated experimental pulse. Finite boundaries reflect outgoing amplitude; the present implementation does not calculate ionization yields or high-harmonic spectra. A 1D soft-core model is not the 3D Coulomb atom.

## References

- J. Crank and P. Nicolson (1947), A practical method for numerical evaluation of solutions of partial differential equations of the heat-conduction type, DOI: 10.1017/S0305004100023197.
- D. J. Griffiths and D. F. Schroeter, Introduction to Quantum Mechanics, 3rd ed., Cambridge University Press (2018).
