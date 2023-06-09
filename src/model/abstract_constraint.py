from abc import ABC, abstractmethod


class AbstractConstraint(ABC):
    @abstractmethod
    def constraint(self, a: list):
        pass

    @abstractmethod
    def dimension(self) -> int:
        pass

    @abstractmethod
    def c_matrix(self, a: list):
        pass
