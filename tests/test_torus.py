import numpy as np
from scipy.linalg import expm
from openatom.torus import TorusGrid, TwoStateTorus, GAP, DIPOLE

def test_torus_geometry():
    g=TorusGrid(8,3,.5)
    x,y,z=g.xyz.T
    assert np.all((np.hypot(x,y)-8)**2+z*z<9)
    assert abs(g.volume/g.exact_volume-1)<.12
    assert (g.kinetic-g.kinetic.T).nnz==0
    assert np.all(g.profile(.3,.1)>0)

def test_field_off_and_norm():
    g=TorusGrid(6,2.5,1)
    p=TwoStateTorus(g,.2)
    state=g.packet([1,1.3,1],.1)
    for _ in range(10):
        state=p.step(state,0)
    np.testing.assert_allclose(p.populations(state),[1,0],atol=1e-12)
    for _ in range(20):
        state=p.step(state,.03*g.profile(.3,.1))
    assert p.populations(state)[1]>0
    np.testing.assert_allclose(p.populations(state).sum(),1,atol=1e-12)

def test_constant_uniform_field_against_exact_internal():
    g=TorusGrid(4,1.9,1)
    p=TwoStateTorus(g,.4)
    state=p.step(g.packet([1,1,1]),.07)
    h=np.array([[0,-DIPOLE*.07],[-DIPOLE*.07,GAP]])
    exact=expm(-1j*.4*h)@np.array([1,0])
    np.testing.assert_allclose(p.populations(state),abs(exact)**2,atol=1e-12)

def test_nonuniform_field_time_refinement_against_full_matrix():
    g=TorusGrid(3,1.4,1,mass=1)
    initial=g.packet([1,1,1])
    field=.2*g.profile(.3,.1)
    n=len(g.xyz)
    H=np.kron(g.kinetic.toarray(),np.eye(2))
    for i,e in enumerate(field):
        H[2*i:2*i+2,2*i:2*i+2]+=np.array([[0,-DIPOLE*e],[-DIPOLE*e,GAP]])
    exact=(expm(-1j*.4*H)@initial.ravel()).reshape(n,2)
    def error(dt):
        p=TwoStateTorus(g,dt)
        state=initial
        for _ in range(round(.4/dt)):
            state=p.step(state,field)
        return np.linalg.norm(state-exact)
    assert error(.05)<error(.1)/3
