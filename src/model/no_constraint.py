from model.abstract_constraint import AbstractConstraint
import numpy as np


class NoConstraint(AbstractConstraint):
    def constraint(self, a: np.ndarray) -> np.ndarray:
        return None

    def dimension(self) -> int:
        return 0

    def c_matrix(self, a: np.ndarray) -> np.ndarray:
        return None
