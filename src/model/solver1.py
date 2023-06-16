import numpy as np
from model.abstract_model import AbstractModel
from model.abstract_solver import AbstractSolver
from model.observation import Observation


class Solver1(AbstractSolver):
    """
    This class is the concrete subclass of AbstractSolver.
    See "Least Square Adjustment with relatively large observation errors,
    inaccurate initial approximations, or both" by Heinrich Eichhorn and Warren
    G. Clay, Mon. Not. R. astr. Soc. (1974) 166, 425-432.
    """
    def __init__(self, m: AbstractModel):
        super().__init__(m)

    def correction(self, a0: np.ndarray, obs: np.ndarray):
        # TODO
        pass

    def covariant(self, a0: np.ndarray, obs: np.ndarray):
        _no: int = len(obs)
        tmp: np.ndarray = np.array([_no + self._nc + self._np,
                                    _no + self._nc + self._np])
        b: np.ndarray = self._model.b_matrix(a0, obs)
        c: np.ndarray = self._model.c_matrix(a0)
        pass

    @staticmethod
    def _matrix_g_inverse(obs: np.ndarray):
        _tmp = []
        for i in range(len(obs)):
            _tmp.append(1 / obs[i].get_sigma() / obs[i].get_sigma())
        return np.diag(_tmp)
