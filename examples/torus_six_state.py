"""Generate six-state pulse CSV and diagnostics; use --matlab-config for matching input."""
import argparse
from openatom.torus_study import PulseConfig, matlab_config, run_pulse, save_study

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",default="torus-study-output")
    p.add_argument("--matlab-config",action="store_true")
    a=p.parse_args()
    config=matlab_config() if a.matlab_config else PulseConfig()
    rows,d=run_pulse(config)
    save_study(a.output,rows,d)
    print(d)
