from dataclasses import replace
import csv
import numpy as np
import pytest
from openatom.torus_study import PulseConfig, run_pulse, save_study
from openatom.benchmarks import compare_pulse_csv

def test_study_csv_and_comparison(tmp_path):
    c=PulseConfig(major=4,minor=1.9,spacing=1,sigma=(1,1,1),window_as=48,step_as=12,center_as=24,sigma_as=12)
    rows,d=run_pulse(c)
    assert len(rows)==5
    assert d["max_norm_error"]<1e-11
    assert rows[-1]["Population_2p0_On"]>0
    assert rows[-1]["Population_2p0_Off"]<1e-20
    save_study(tmp_path/"a",rows,d)
    a=tmp_path/"a"/"pulse_full_observables.csv"
    # Synthetic perturbation only: tests comparison mechanics, not MATLAB agreement.
    modified=[dict(r) for r in rows]
    modified[-1]["Population_2p0_On"]+=.001
    save_study(tmp_path/"b",modified,d)
    b=tmp_path/"b"/"pulse_full_observables.csv"
    errors=compare_pulse_csv(a,b)
    np.testing.assert_allclose(errors["Population_2p0_On"]["max_abs"],.001)
    assert errors["Norm_On"]["max_abs"]==0
    modified[-1]["Time_as"]+=1
    save_study(tmp_path/"b",modified,d)
    with pytest.raises(ValueError):
        compare_pulse_csv(a,b)
