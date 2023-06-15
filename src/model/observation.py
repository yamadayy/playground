import numpy as np


class Observation:
    def __init__(self, v: float, s: float, p: np.ndarray):
        self.__value: float = v
        self.__sigma: float = s
        self.__params = None
        if not (p is None):
            self.__params = p.copy()

    def get_value(self) -> float:
        return self.__value

    def get_sigma(self) -> float:
        return self.__sigma

    def get_params(self) -> np.ndarray:
        return self.__params

    def set_params(self, p: np.ndarray):
        self.__params = p.copy()
