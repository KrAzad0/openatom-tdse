import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm
from openatom.multilevel import NL, radial, hydrogen_basis, MultiStateTorus
from openatom.torus import TorusGrid, DIPOLE, GAP

def test_basis_analytic():
    e,d,r=hydrogen_basis()
    np.testing.assert_allclose(e[2],GAP,atol=1e-13)
    np.testing.assert_allclose(d[0,2],DIPOLE,atol=1e-10)
    np.testing.assert_allclose(d,d.T,atol=1e-13)
    np.testing.assert_allclose(r,r.T,atol=1e-13)
    for i,(n,l) in enumerate(NL):
        np.testing.assert_allclose(quad(lambda x: radial(n,l,x)**2*x*x,0,np.inf)[0],1,atol=1e-9)
        for j,(_,ll) in enumerate(NL):
            if abs(l-ll)!=1:
                assert d[i,j]==0

def test_field_off_and_observables():
    g=TorusGrid(4,1.9,1)
    e,d,r=hydrogen_basis()
    p=MultiStateTorus(g,.2,e,d)
    v=p.initial([1,1,1])
    for _ in range(3):
        v=p.step(v,0)
    obs=p.observe(v,r)
    np.testing.assert_allclose(obs["populations"],[1,0,0,0,0,0],atol=1e-12)
    np.testing.assert_allclose(obs["relative_z"],0,atol=1e-12)
    v=p.step(v,.1*g.profile(.2))
    obs=p.observe(v,r)
    np.testing.assert_allclose(obs["norm"],1,atol=1e-12)
    np.testing.assert_allclose(obs["electron_z"]-obs["proton_z"],obs["relative_z"],atol=1e-12)
    assert obs["populations"][2]>0

def test_uniform_field_temporal_order():
    g=TorusGrid(4,1.9,1)
    e,d,_=hydrogen_basis()
    h=np.diag(e)-.1*d
    initial=np.eye(6)[0]
    exact=abs(expm(-1j*h)@initial)**2
    def error(dt):
        p=MultiStateTorus(g,dt,e,d)
        v=p.initial([1,1,1])
        for _ in range(round(1/dt)):
            v=p.step(v,.1)
        return np.linalg.norm(p.observe(v)["populations"]-exact)
    assert error(.1)<error(.2)/3
