import numpy as np

from jasmine_toolkit.optimize.abstract_model import AbstractModel
from jasmine_toolkit.optimize.abstract_solver import AbstractSolver
from jasmine_toolkit.datamodel.toy.free_fall_model import FreeFallModel
from jasmine_toolkit.optimize.least_square_fit import LeastSquareFit
from jasmine_toolkit.optimize.observation import Observation
from jasmine_toolkit.optimize.solver1 import Solver1


def test_free_fall():
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
    for i in range(n):
        print(obs[i].get_value())
    solver: AbstractSolver = Solver1(model)
    lsf: LeastSquareFit = LeastSquareFit(model, solver)
    for i in range(n):
        lsf.add_observation(obs[i])
    lsf.solve(a)
