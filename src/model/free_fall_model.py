import numpy as np

from model.abstract_model import AbstractModel
from model.abstract_solver import AbstractSolver
from model.least_square_fit import LeastSquareFit
from model.observation import Observation
from model.solver1 import Solver1


class FreeFallModel(AbstractModel):

    def __init__(self):
        super().__init__(3)

    def model(self, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        return a[0] * b[0] * b[0] + a[1] * b[0] + a[2]

    def b_matrix(self, a: np.ndarray, o: np.ndarray) -> np.ndarray:
        tmp = np.zeros((len(o), 3))
        for i in range(len(o)):
            tmp[i][0] = o[i] * o[i]
            tmp[i][1] = o[i]
            tmp[i][2] = 1.
        return tmp


if __name__ == '__main__':
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
    """
    model: AbstractModel = FreeFallModel()
    print('number of parameters:  ' + str(model.num_param()))
    print('number of constraints: ' + str(model.num_constraint()))
    n = 11
    oo = []
    aa = [1.0, 0.0, 0.0]
    error = 0.01
    r = np.random.normal(loc=0, scale=error, size=n)
    for ii in range(n):
        t = 0.1 * ii
        oo.append(model.model(aa, [t]) + 0 * r[ii])
    print('observed values are ' + str(oo))
    mat_b = model.b_matrix(aa, oo)
    print('B = ' + str(mat_b))
    print('C = ' + str(model.c_matrix(aa)))
    a_try = [0.8, 0.0, 0.0]
    sigma = np.identity(n) * error * error
    F0 = np.zeros((n, 1))
    for j in range(10):
        for ii in range(n):
            F0[ii][0] = model.model(a_try, [0.1 * ii]) - oo[ii]
        sigma_param = np.linalg.inv(np.dot(mat_b.T, np.dot(np.linalg.inv(sigma),
                                                           mat_b)))
        delta = - np.dot(sigma_param, np.dot(mat_b.T,
                                             np.dot(np.linalg.inv(sigma), F0)))
        for ii in range(3):
            a_try[ii] = a_try[ii] + delta[ii][0]
        print(str(delta) + " : " + str(a_try))
    """
