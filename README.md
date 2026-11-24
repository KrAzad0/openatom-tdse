# OpenAtom-TDSE

An early-stage open-source Python toolkit for computational atomic physics and time-dependent quantum mechanics. The aim is reproducible numerical examples for students and researchers working on light–matter interactions.

## Current scope

- Radial hydrogen-like Coulomb eigenstates via a tridiagonal finite-difference Hamiltonian.
- One-dimensional real-potential midpoint Crank–Nicolson propagation.
- Illustrative Gaussian-driven soft-core atom example.
- Analytic spectrum, orthogonality, convergence and discrete unitarity tests.

This is development software, not an experimentally validated research package. Three-dimensional propagation, attosecond observables and toroidal confinement are planned, not implemented.

## Install and run

```sh
git clone https://github.com/KrAzad0/openatom-tdse.git
cd openatom-tdse
python -m pip install -e ".[dev]"
python -m pytest
python examples/hydrogen_spectrum.py
python examples/softcore_pulse.py
```

See [numerical theory and limitations](docs/theory.md), [contributing](CONTRIBUTING.md) and [roadmap](ROADMAP.md).

## License and citation

MIT; see LICENSE. Cite the repository and exact commit used (see CITATION.cff). No published package or external adoption is claimed.
