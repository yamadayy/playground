import math

import numpy as np

from jasmine_toolkit.optimize.abstract_model import AbstractModel
from jasmine_toolkit.optimize.abstract_solver import AbstractSolver
from jasmine_toolkit.optimize.observation import Observation


class LeastSquareFit:
    """
    See "Least Square Adjustment with relatively large observation errors,
    inaccurate initial approximations, or both" by Heinrich Eichhorn and Warren
    G. Clay, Mon. Not. R. astr. Soc. (1974) 166, 425-432.
    """
    def __init__(self, m: AbstractModel, s: AbstractSolver):
        self.__model: AbstractModel = m
        self.__solver: AbstractSolver = s
        self.__observations: np.ndarray = np.empty(0)

    def dimension(self):
        """

        :return: number of parameters which we would like to estimate.
        """
        return self.__model.num_param()

    def add_observation(self, o: Observation):
        """

        :param o: A single observation (An Observation object)
        :return: none
        """
        self.__observations = np.append(self.__observations, o)

    def get_model(self) -> AbstractModel:
        """

        :return: model which we use.
        """
        return self.__model

    def covariant(self, a0: np.ndarray) -> np.ndarray:
        """
        This routine simply call the method with the same name in solver class.
        :param a0: trial values of parameters which we would like to estimate
        :return: covariant matrix of estimate parameters, Sigam.
        """
        return self.__solver.covariant(a0, self.__observations)

    def solve(self, a0: np.ndarray, threshold: float = 0.01,
              max_iter: int = 10) -> np.ndarray:
        """

        :param a0: trial values of parameters which we would like to estimate
        :param threshold: threshold of the residual where the iteration stops.
        :param max_iter: maximum number of iteration
        :return: corrected value of {a}.
        """
        tmp = np.zeros(len(a0))
        cor = self.__solver.correction(a0, self.__observations)
        tmp = a0 + cor
        return tmp

    def get_S(self, a0: np.ndarray) -> float:
        """

        :param a0: trial values of parameters which we would like to estimate
        :return: residual.
        """
        tmp = 0.0
        for i in range(len(self.__observasions)):
            o: Observation = self.__observasions.get(i)
            tmp = tmp + math.pow(o.get_value()
                                 - self.__model.model(a0, o.get_params()), 2.0)
        return tmp
