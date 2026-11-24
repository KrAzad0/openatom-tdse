"""Small field-on/off torus example; illustrative, not a MATLAB reproduction."""
import numpy as np
from openatom import gaussian_pulse
from openatom.torus import TorusGrid, TwoStateTorus, ATOMIC_TIME_AS, GAP

grid=TorusGrid(26,8,4)
dt=12/ATOMIC_TIME_AS
p=TwoStateTorus(grid,dt)
on=off=grid.packet([2.4,4,2.4],.1)
profile=grid.profile(.25)
print("Interior nodes:",len(grid.xyz))
print("Relative grid volume error:",grid.volume/grid.exact_volume-1)
for k in range(120):
    field=gaussian_pulse((k+.5)*dt,.015,GAP,650/ATOMIC_TIME_AS,180/ATOMIC_TIME_AS)
    on=p.step(on,field*profile)
    off=p.step(off,0)
print("Field-on populations:",p.populations(on))
print("Field-off populations:",p.populations(off))
print("Field-on norm error:",abs(p.populations(on).sum()-1))
