import math

from jasmine_toolkit.optimize.abstract_model import AbstractModel
from jasmine_toolkit.optimize.observation import Observation
from jasmine_toolkit.optimize.solver1 import Solver1
from itertools import product

from jasmine_toolkit.datamodel.toy.free_fall_model import FreeFallModel
import numpy as np


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


def test_large_matrix():
    model: AbstractModel = FreeFallModel()
    a = np.array([1.0, 0.0, 0.0])
    obs = np.empty(0)
    n = 11
    error = 0.01
    r = np.random.normal(loc=0, scale=error, size=n)
    for i in range(n):
        b = np.array([i * 0.1])
        o = Observation(model.model(a, b) + r[i], error, b)
        obs = np.append(obs, np.array([o]))
    solver: Solver1 = Solver1(model)
    tmp = solver._large_matrix(a, obs)
    for i in range(n):
        assert math.isclose(tmp[i][n], i * i * 0.01, abs_tol=1.0e-7)
        assert math.isclose(tmp[i][n + 1], i * 0.1, abs_tol=1.0e-7)
        assert math.isclose(tmp[i][n + 2], 1, abs_tol=1.0e-7)


def test_correction():
    model: AbstractModel = FreeFallModel()
    a = np.array([1.0, 0.0, 0.0])
    obs = np.empty(0)
    n = 11
    error = 0.001
    r = np.random.normal(loc=0, scale=error, size=n)
    for i in range(n):
        b = np.array([i * 0.1])
        o = Observation(model.model(a, b) + r[i], error, b)
        obs = np.append(obs, np.array([o]))
    solver: Solver1 = Solver1(model)
    tmp = solver.correction([0., 0., 0.], obs)
    assert math.isclose(tmp[0], 1.0, abs_tol=0.05)
    assert math.isclose(tmp[1], 0.0, abs_tol=0.05)
    assert math.isclose(tmp[2], 0.0, abs_tol=0.05)
