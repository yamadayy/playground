from model.observation import Observation
from model.solver1 import Solver1


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
    for i in range(3):
        for j in range(3):
            if not i == j:
                assert g[i][j] == 0
