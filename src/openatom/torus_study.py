"""Reproducible pulse studies for the six-state composite atom."""
from dataclasses import dataclass, asdict
import csv
import json
from pathlib import Path
import numpy as np
from .torus import TorusGrid, ATOMIC_TIME_AS, GAP
from .tdse import gaussian_pulse
from .multilevel import LABELS, MultiStateTorus, hydrogen_basis

@dataclass(frozen=True)
class PulseConfig:
    major: float=12.
    minor: float=5.
    spacing: float=2.
    sigma: tuple=(2.,3.,2.)
    momentum: float=.07
    amplitude: float=.022
    angular: float=.35
    vertical: float=.10
    center_as: float=900.
    sigma_as: float=220.
    window_as: float=1800.
    step_as: float=24.
    carrier_multiplier: float=1.
    phase: float=0.

def matlab_config():
    """Original six-state MATLAB pulse parameters (costly 3D grid)."""
    return PulseConfig(major=85,minor=32,spacing=6,sigma=(10,14,10),step_as=12)

def run_pulse(config=PulseConfig()):
    """Return numerical time series and diagnostics; run on/off with same discretization."""
    if not np.isfinite([config.window_as,config.step_as,config.sigma_as,config.amplitude,config.center_as,config.phase,config.carrier_multiplier]).all() or min(config.window_as,config.step_as,config.sigma_as)<=0 or config.amplitude<0:
        raise ValueError("invalid pulse parameters")
    grid=TorusGrid(config.major,config.minor,config.spacing)
    times=np.linspace(0,config.window_as,int(np.ceil(config.window_as/config.step_as))+1)
    dt=(times[1]-times[0])/ATOMIC_TIME_AS
    energies,d,radius=hydrogen_basis()
    propagator=MultiStateTorus(grid,dt,energies,d)
    on=off=propagator.initial(config.sigma,config.momentum)
    profile=grid.profile(config.angular,config.vertical)
    def field(t):
        return gaussian_pulse(t/ATOMIC_TIME_AS,config.amplitude,config.carrier_multiplier*GAP,config.center_as/ATOMIC_TIME_AS,config.sigma_as/ATOMIC_TIME_AS,config.phase)
    rows=[]
    for i,t in enumerate(times):
        if i:
            on=propagator.step(on,field((times[i-1]+t)/2)*profile)
            off=propagator.step(off,0)
        row={"Time_as":float(t),"Eref_au":float((1+config.angular)*field(t))}
        for suffix,state in (("On",on),("Off",off)):
            obs=propagator.observe(state,radius)
            row.update({f"Population_{label}_{suffix}":float(pop) for label,pop in zip(LABELS,obs["populations"])})
            row.update({f"Norm_{suffix}":obs["norm"],f"RelativeZ_{suffix}_a0":obs["relative_z"],
                        f"MeanRadius_{suffix}_a0":obs["mean_radius"],f"Xcm_{suffix}_a0":float(obs["com"][0]),
                        f"Ycm_{suffix}_a0":float(obs["com"][1]),f"Zcm_{suffix}_a0":float(obs["com"][2])})
            if suffix=="On":
                row["ElectronZ_On_a0"]=obs["electron_z"]
                row["ProtonZ_On_a0"]=obs["proton_z"]
        rows.append(row)
    diagnostics={"config":asdict(config),"nodes":len(grid.xyz),"actual_step_as":float(times[1]-times[0]),
                 "relative_grid_volume_error":grid.volume/grid.exact_volume-1,
                 "max_norm_error":max(abs(r[f"Norm_{s}"]-1) for r in rows for s in ("On","Off")),
                 "final_populations":{k:rows[-1][f"Population_{k}_On"] for k in LABELS}}
    return rows,diagnostics

def save_study(directory,rows,diagnostics):
    folder=Path(directory)
    folder.mkdir(parents=True,exist_ok=True)
    with (folder/"pulse_full_observables.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (folder/"diagnostics.json").write_text(json.dumps(diagnostics,indent=2)+"\n")
