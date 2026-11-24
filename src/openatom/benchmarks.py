"""CSV comparison utilities. No external reference results are bundled."""
import csv
import numpy as np
from .multilevel import LABELS

def read_csv(path):
    with open(path,newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

def compare_pulse_csv(reference,candidate):
    """Compare matched time grids; return max absolute and RMS errors per column.

    Both simulations must use the same physical parameters, grid, basis phases
    and drive sign. No interpolation or tolerance-based pass is implicit.
    """
    a,b=read_csv(reference),read_csv(candidate)
    if not a or not b or len(a)!=len(b):
        raise ValueError("nonempty CSVs with equal sample counts required")
    def column(rows,key):
        try:
            v=np.array([float(r[key]) for r in rows])
        except (KeyError,TypeError,ValueError) as exc:
            raise ValueError(f"invalid or missing column {key}") from exc
        if not np.isfinite(v).all():
            raise ValueError(f"nonfinite column {key}")
        return v
    ta,tb=column(a,"Time_as"),column(b,"Time_as")
    if np.any(np.diff(ta)<=0) or not np.allclose(ta,tb,rtol=0,atol=1e-7):
        raise ValueError("time grids must match and increase")
    keys=["Eref_au"]+[f"Population_{k}_{s}" for k in LABELS for s in ("On","Off")]+[f"Norm_{s}" for s in ("On","Off")]
    keys += [k for k in ("RelativeZ_On_a0","RelativeZ_Off_a0","MeanRadius_On_a0","MeanRadius_Off_a0","ElectronZ_On_a0","ProtonZ_On_a0") if k in a[0] and k in b[0]]
    out={}
    for k in keys:
        delta=column(a,k)-column(b,k)
        out[k]={"max_abs":float(np.max(abs(delta))),"rms":float(np.sqrt(np.mean(delta**2)))}
    return out
