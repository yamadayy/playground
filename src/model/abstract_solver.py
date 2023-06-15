from abc import ABC, abstractmethod
from model.abstract_model import AbstractModel
import numpy as np


class AbstractSolver(ABC):
    def __init__(self, m: AbstractModel):
        self.model: AbstractModel = m
        self.n = self.model.num_param()

    @abstractmethod
    def correlation(self, a: np.ndarray, da: np.ndarray, o: np.ndarray):
        pass

    @abstractmethod
    def covariant(self, a: np.ndarray, da: np.ndarray, o: np.ndarray):
        pass
