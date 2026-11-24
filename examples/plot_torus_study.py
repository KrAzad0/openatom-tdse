"""Plot saved six-state pulse populations and field; works headless."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openatom.benchmarks import read_csv
from openatom.multilevel import LABELS

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("csv")
    p.add_argument("--output",default="torus-plots")
    a=p.parse_args()
    rows=read_csv(a.csv)
    t=[float(r["Time_as"]) for r in rows]
    out=Path(a.output)
    out.mkdir(parents=True,exist_ok=True)
    fig,ax=plt.subplots()
    for label in LABELS:
        ax.plot(t,[float(r[f"Population_{label}_On"]) for r in rows],label=label)
    ax.set(xlabel="Time (as)",ylabel="State population",title="Six-state toroidal hydrogen")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out/"populations.png",dpi=160)
    plt.close(fig)
    fig,ax=plt.subplots()
    ax.plot(t,[float(r["Eref_au"]) for r in rows])
    ax.set(xlabel="Time (as)",ylabel="Reference drive (a.u.)")
    fig.tight_layout()
    fig.savefig(out/"drive.png",dpi=160)
    plt.close(fig)
