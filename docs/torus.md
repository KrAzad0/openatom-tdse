# Experimental toroidal composite-atom solver

Inspired by the maintainer's November 2026 two-state and six-state MATLAB research models. This is an independent Python implementation of the two-state method, not a port of every MATLAB feature or a reproduction of its numerical results. No original MATLAB files or unpublished data are included.

The solid torus is defined by (sqrt(X^2+Y^2)-R)^2+Z^2<a^2, with volume 2*pi^2*R*a^2. The Cartesian seven-point kinetic stencil retains the full diagonal and removes exterior neighbors (zero exterior wavefunction). The mass is the electron plus proton mass in atomic units; the internal gap is 3*mu/8, and the magnitude of the 1s–2p_z position matrix element is 256/(243*sqrt(2)*mu).

Expand the composite state in 1s and 2p_z. The local internal Hamiltonian has diagonals 0 and gap and off-diagonals -d*E, matching the original signed drive convention. The physical electric dipole of hydrogen is minus its electron-proton relative coordinate; mapping the drive sign to a laboratory electric field requires a consistent basis phase and dipole convention. Populations are invariant under reversal of the off-diagonal sign for this initial state.

The step uses kinetic half-step Cayley transforms and the exact local 2x2 exponential at the temporal midpoint. The discrete norm is h^3 times the sum of squared amplitudes. The tests check symmetry, geometry, field-off populations, norm, exact uniform-field internal populations, and time refinement against a full small-matrix exponential for a spatially varying static field.

Run: python examples/torus_two_state.py

Limitations: a staircase approximation to the torus wall; spatial grid and time convergence must be assessed for actual research parameters. Walls act on COM only, not separate particles. Free-space internal states, no ionization, no self-consistent Maxwell solution, no continuum or six-state dynamics yet. The example is deliberately coarse for demonstration. Norm conservation is not proof of convergence. Do not infer experimental predictions from this example.

Next: systematic spatial and temporal convergence, six-state hydrogen dipoles and comparison with independently rerun MATLAB benchmarks.
