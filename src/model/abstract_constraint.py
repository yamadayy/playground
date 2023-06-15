from abc import ABC, abstractmethod
import numpy as np


class AbstractConstraint(ABC):
    @abstractmethod
    def constraint(self, a: np.ndarray) -> np.ndarray:
        pass

    @abstractmethod
    def dimension(self) -> int:
        pass

    @abstractmethod
    def c_matrix(self, a: np.ndarray) -> np.ndarray:
        pass
