from abc import ABC, abstractmethod
from model.no_constraint import NoConstraint
import numpy as np


class AbstractModel(ABC):
    def __init__(self, n: int):
        self.__n = n
        self.__constraint = NoConstraint()

    @abstractmethod
    def model(self, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        pass

    @abstractmethod
    def b_matrix(self, a: list, o: list):
        pass

    def num_param(self):
        return self.__n

    def num_constraint(self):
        return self.__constraint.dimension()

    def c_matrix(self, a: list):
        return self.__constraint.c_matrix(a)

    def constraint(self, a: list):
        return self.__constraint.constraint(a)
