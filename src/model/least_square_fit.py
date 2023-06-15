import numpy as np

from model.abstract_model import AbstractModel
from model.abstract_solver import AbstractSolver
from model.observation import Observation


class LeastSquareFit:
    def __init__(self, m: AbstractModel, s: AbstractSolver):
        self.__model = m
        self.__solver = s
        self.__observations: np.ndarray = None

    def dimension(self):
        return self.__model.num_param()

    def add_observation(self, o: Observation):
        self.__observations = np.append(self.__observations, o)

    def get_model(self) -> AbstractModel:
        return self.__model

    def covariant(self, a: np.ndarray, da: np.ndarray) -> np.ndarray:
        return self.__solver.covariant(a, da, self.__observations)

    def solve(self, a: np.ndarray, da: np.ndarray) -> np.ndarray:
        # TODO implements this function.
        pass

    def get_S(self, a: np.ndarray) -> float:
        # TODO implements this function
        return 0.0
