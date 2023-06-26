import numpy as np

from model.abstract_model import AbstractModel
from model.abstract_solver import AbstractSolver


class Solver1(AbstractSolver):
    """
    This class is a concrete subclass of AbstractSolver.
    See "Least Square Adjustment with relatively large observation errors,
    inaccurate initial approximations, or both" by Heinrich Eichhorn and Warren
    G. Clay, Mon. Not. R. astr. Soc. (1974) 166, 425-432.
    """

    def __init__(self, m: AbstractModel):
        super().__init__(m)

    def correction(self, a0: np.ndarray, obs: np.ndarray):
        # TODO (under construction)
        large_matrix = self._large_matrix(a0, obs)
        large_matrix_inverse = - np.linalg.inv(large_matrix)
        g0 = np.empty(0)
        for i in range(len(obs)):
            g0 = np.append(g0, self._model.model(a0, obs[i].get_params()))
        tmp_f = np.zeros(len(obs) + self._nc + self._np)
        _no = len(obs)
        tmp_f[0:_no] = g0[0:_no]
        for i in range(_no):
            tmp_f[i] = tmp_f[i] + obs[i].get_value()
        if self._nc > 0:
            h0 = self._model.constraint(a0)
            tmp_f[_no:_no + self._nc] = h0[0:self._nc]
        ans = large_matrix_inverse.dot(tmp_f)
        return ans[_no + self._nc:_no + self._nc + self._np]

    def covariant(self, a0: np.ndarray, obs: np.ndarray):
        """
        Covariant of the estimate parameter appears in the right bottom sub
        matrix of the inverse of large matrix.
        :param a0: trial values of parameters which we would like to estimate
        :param obs: a set of observations
        :return: covariant of estimate parameters, Sigma
        """
        # TODO need optimize
        tmp = - np.lilnalg.inv(self._large_matrix(a0, obs))
        _no: int = len(obs)
        return tmp[_no + self._nc:_no + self._nc + self._np,
                   _no + self._nc:_no + self._nc + self._np]

    @staticmethod
    def _matrix_g_inverse(obs: np.ndarray):
        """
        A matrix sigma is covariant matrix of observation.
        G = X sigma X^T, and the problem is formulated as a matrix X be unity.
        So, this routine returns sigma^(-1).
        :param obs: a set of observations
        :return: inverse of matrix G
        """
        _tmp = []
        for i in range(len(obs)):
            _tmp.append(1 / obs[i].get_sigma() / obs[i].get_sigma())
        return np.diag(_tmp)

    def _large_matrix(self, a0: np.ndarray, obs: np.ndarray) -> np.ndarray:
        """
        In the paper, the Least Square Problem is written in matrix form
         |X sigma X^T 0   B | |Lambda_0|   |G_0|   |0|
         |     0      0   C | |Lambda_1| + |H_0| = |0|
         | B^T       C^T  0 | | Delta  |   | 0 |   |0|.
         This routine returns the coefficient matrix.

        :param a0: trial values of parameters which we would like to estimate
        :param obs: a set of observations
        :return: matrix form of the LSF problem.
        """
        _no: int = len(obs)
        tmp: np.ndarray = np.zeros([_no + self._nc + self._np,
                                    _no + self._nc + self._np])
        b: np.ndarray = self._model.b_matrix(a0, obs)
        c: np.ndarray = self._model.c_matrix(a0)
        gd: np.ndarray = np.empty(0)
        for i in range(len(obs)):
            gd = np.append(gd, obs[i].get_sigma() * obs[i].get_sigma())
        g = np.diag(gd)
        tmp[0:_no, 0:_no] = g[0:_no, 0:_no]
        tmp[0:_no, _no + self._nc:_no + self._nc + self._np] \
            = b[0:_no, 0:self._np]
        tmp[_no + self._nc:_no + self._nc + self._np, 0:_no] \
            = b[0:_no, 0:self._np].T
        if self._nc > 0:
            tmp[_no:_no + self._nc, _no + self._nc:_no + self._nc + self._np] \
                = c[0:self._nc, 0:self._np]
            tmp[_no + self._nc:_no + self._nc + self._np, _no:_no + self._nc] \
                = c[0:self._nc, 0:self._np].T
        return tmp
