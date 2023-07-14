from jasmine_toolkit.optimize.observation import Observation
from jasmine_toolkit.optimize.solver1 import Solver1
from itertools import product


def test_g_inverse():
    tmp = []
    param = [1.0]
    tmp.append(Observation(1.0, 1.0, param))
    tmp.append(Observation(1.0, 0.5, param))
    tmp.append(Observation(1.0, 0.25, param))
    g = Solver1._matrix_g_inverse(tmp)
    assert g[0][0] == 1
    assert g[1][1] == 4
    assert g[2][2] == 16
    for i, j in product(range(3), range(3)):
        assert g[i][j] == 0 or i == j
