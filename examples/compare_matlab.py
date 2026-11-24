"""Compare two pulse CSV files on identical time samples."""
import argparse
import json
from openatom.benchmarks import compare_pulse_csv

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("matlab_csv")
    p.add_argument("python_csv")
    a=p.parse_args()
    print(json.dumps(compare_pulse_csv(a.matlab_csv,a.python_csv),indent=2))
