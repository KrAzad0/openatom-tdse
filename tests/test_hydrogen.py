import numpy as np
import pytest
from openatom import hydrogen_radial

@pytest.mark.parametrize("l", [0, 1, 2])
def test_spectrum_and_orthogonality(l):
    r, e, u = hydrogen_radial(l=l, r_max=100, points=2500, states=3)
    n = np.arange(l+1, l+4)
    np.testing.assert_allclose(e, -0.5/n**2, atol=5e-4, rtol=0)
    np.testing.assert_allclose((r[1]-r[0]) * u.T @ u, np.eye(3), atol=1e-12)

def test_ground_state_convergence():
    _, coarse, _ = hydrogen_radial(points=400, states=1)
    _, fine, _ = hydrogen_radial(points=800, states=1)
    assert abs(fine[0]+0.5) < abs(coarse[0]+0.5)/3

def test_invalid_l():
    with pytest.raises(ValueError):
        hydrogen_radial(l=-1)
