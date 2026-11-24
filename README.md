# OpenAtom-TDSE

An early-stage MIT-licensed Python toolkit for computational atomic physics and time-dependent quantum mechanics. The aim is reproducible numerical examples for students and researchers.

## Current scope

- Radial hydrogen-like Coulomb eigenstates via tridiagonal finite differences.
- One-dimensional real-potential midpoint Crank–Nicolson propagation and a Gaussian-driven soft-core atom example.
- Experimental solid-torus center-of-mass dynamics with two or six internal bound hydrogen states.
- Field-on/off pulse studies, CSV observables, parameter sweeps, numerical refinement tools and plotting.
- Analytical spectrum, orthogonality, dipole selection, norm and finite-matrix reference tests.

This is development software, not an experimentally validated research package. The torus wall confines the atomic center of mass only; the field is prescribed and the internal basis excludes ionization. Research-scale convergence and independent MATLAB agreement remain to be established.

## Install and run

```sh
git clone https://github.com/KrAzad0/openatom-tdse.git
cd openatom-tdse
python -m pip install -e ".[dev]"
python -m pytest
python examples/hydrogen_spectrum.py
python examples/softcore_pulse.py
python examples/torus_two_state.py
python examples/torus_six_state.py
python examples/plot_torus_study.py torus-study-output/pulse_full_observables.csv
python examples/torus_parameter_study.py
```

See [numerical theory](docs/theory.md), [two-state torus model](docs/torus.md), [six-state and benchmark protocol](docs/torus_benchmarks.md), [software-paper draft](docs/software_paper_draft.md), [contributing](CONTRIBUTING.md) and [roadmap](ROADMAP.md).

## License and citation

MIT; see LICENSE. Cite the repository and exact commit used (see CITATION.cff). The initial development release is v0.1.0-dev.0; subsequent development may include additional features not in that release. No PyPI publication, independent adoption or archival DOI is claimed.
