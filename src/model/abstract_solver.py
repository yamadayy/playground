from abc import ABC, abstractmethod
from model.abstract_model import AbstractModel


class AbstractSolver(ABC):
    def __init__(self, m: AbstractModel):
        self.model: AbstractModel = m
        self.n = self.model.num_param()

    @abstractmethod
    def correlation(self, a: list, da: list, o: list):
        pass

    @abstractmethod
    def covariant(self, a: list, da: list, o: list):
        pass
