"""Small deterministic demonstration of drive-amplitude and temporal/spatial refinement."""
from dataclasses import replace
import argparse
import json
from pathlib import Path
from openatom.torus_study import PulseConfig, run_pulse

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",default="torus-parameter-study.json")
    a=p.parse_args()
    base=PulseConfig(major=8,minor=3,spacing=1.5,sigma=(1.5,2,1.5),step_as=48)
    def case(**kwargs):
        _,d=run_pulse(replace(base,**kwargs))
        return d
    result={"amplitude_sweep":[case(amplitude=x) for x in (0,.01,.02,.04)],
            "time_refinement":[case(step_as=x) for x in (48,24,12)],
            "space_refinement":[case(spacing=x) for x in (2,1.5,1)]}
    for key in ("time_refinement","space_refinement"):
        seq=result[key]
        result[key+"_successive_max_population_differences"]=[max(abs(a["final_populations"][s]-b["final_populations"][s]) for s in a["final_populations"]) for a,b in zip(seq,seq[1:])]
    Path(a.output).write_text(json.dumps(result,indent=2)+"\n")
    print("Wrote",a.output)
