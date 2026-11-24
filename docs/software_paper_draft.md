# OpenAtom-TDSE: reproducible numerical models for atomic and toroidal quantum dynamics

**Working software-paper draft — not submitted, peer reviewed or assigned a DOI.**

## Summary

OpenAtom-TDSE is a Python toolkit for radial Coulomb eigenstates, one-dimensional time propagation, and experimental composite-atom motion inside a solid toroidal domain. The torus implementation expands the internal hydrogen state in a finite bound-state basis and evolves the center-of-mass amplitudes with sparse finite differences and symmetric time splitting.

## Statement of need

Students exploring ultrafast atomic dynamics benefit from explicit numerical approximations, reproducible examples and tests against analytical reference results. OpenAtom-TDSE aims to expose these approximations and distinguish numerical accuracy from physical model validity. Independent adoption and comparison with existing packages remain to be documented.

## Methods and verification

The package uses a tridiagonal radial Coulomb Hamiltonian, midpoint Crank–Nicolson propagation, a seven-point zero-exterior COM stencil, and two- or six-state internal dynamics. Tests target analytic Coulomb energies, orthogonality, dipole selection rules, discrete norm conservation and temporal refinement against exact finite-matrix evolution. Cross-language MATLAB comparison tooling, parameter sweeps and visualizations are included; numerical agreement with the prior MATLAB runs is not yet established.

## Reproducibility and limitations

See torus_benchmarks.md for commands, parameter matching and convergence protocol. The finite internal basis excludes continuum and ionization; the field is prescribed and toroidal walls constrain only COM. Additional independent verification is required before research predictions are asserted.

## References

- Crank, J. and Nicolson, P. (1947), A practical method for numerical evaluation of solutions of partial differential equations of the heat-conduction type, doi:10.1017/S0305004100023197.
- Griffiths, D. J. and Schroeter, D. F. (2018), Introduction to Quantum Mechanics, 3rd edition, Cambridge University Press.

## Archiving checklist

Review code and manuscript; complete actual numerical comparisons, report measured results and software versions; create a versioned release; connect the repository to Zenodo and archive that release; insert the real DOI in citation metadata and the manuscript. Do not claim a DOI until issued. Verify authorship and licensing before submission.
