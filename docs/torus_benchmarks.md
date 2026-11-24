# Six-state validation and reproducibility protocol

The Python implementation uses the same six free-space m=0 hydrogen basis states, reduced mass, seven-point toroidal COM stencil and nested symmetric splitting as the maintainer's advanced MATLAB model. The signed drive convention H=H0-E*D is retained for comparison. The physical laboratory-field sign and orbital phase must be consistent before interpreting signed dipoles.

Install with python -m pip install -e ".[dev]"; run python -m pytest.

## Cross-language comparison (contribution 4)

1. Rerun the original MATLAB six-state code and retain pulse_full_observables.csv plus exact parameters, MATLAB version and diagnostic log locally.
2. Run python examples/torus_six_state.py --matlab-config --output matched-python (large grid; considerable memory and runtime).
3. Run python examples/compare_matlab.py path/to/MATLAB/pulse_full_observables.csv matched-python/pulse_full_observables.csv.
4. Verify grid counts, parameters, units and basis phases. The script compares matching times only and reports errors without inventing a pass threshold. Establish tolerances from refinement before interpreting differences.

No independent MATLAB reference data have been published or used to claim numerical agreement. Only the pulse window is implemented here, not the MATLAB long-duration transport or density plots. The two-state module remains a separate approximation.

## Pulse study and convergence (contributions 5–6)

python examples/torus_parameter_study.py --output study.json

The script records final populations, norm errors, grid volumes and successive population differences across drive amplitudes, time steps and spatial spacings. These are coarse demonstration settings; no claim of physical convergence follows. Refine time at fixed grid, space at fixed small time step, and enlarge basis/domain as appropriate. Staircase boundaries can reduce spatial convergence order. For publication also compare COM observables, wall probability, energy and full basis truncation effects.

## Visualization (contribution 7)

python examples/torus_six_state.py --output torus-study-output
python examples/plot_torus_study.py torus-study-output/pulse_full_observables.csv

Save the exact commit, dependency versions, parameters and diagnostics with each run. The plots are generated from actual computed CSV data, not illustrative values. Default output paths are local; do not commit unpublished MATLAB data.

## Scientific limitations

No ionization continuum, higher m states, spin, radiation pressure, self-consistent Maxwell field or individual electron/proton wall boundary. The wall acts on COM only. The Gaussian field is transient and is spatially imposed; the model does not establish physical feasibility of such a field in a nanoscale torus. Norm preservation is necessary but insufficient.
