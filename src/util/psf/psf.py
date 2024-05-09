from abc import ABC, ABCMeta, abstractmethod


class Psf(metaclass=ABCMeta):
    @abstractmethod
    def get(self, x: float, y: float):
        pass
