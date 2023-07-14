import numpy as np


class Observation:
    """
    This class is data holder of single observation.
    """
    def __init__(self, v: float, s: float, p: np.ndarray):
        """

        :param v: The value which we get by observation
        :param s: Expected standard deviation of observation error
        :param p: additional parameters like time etc.
        """
        self.__value: float = v
        self.__sigma: float = s
        self.__params = None
        if not (p is None):
            self.__params = p.copy()

    def get_value(self) -> float:
        """

        :return: observed value
        """
        return self.__value

    def get_sigma(self) -> float:
        """

        :return: expected standard deviation of observation error
        """
        return self.__sigma

    def get_params(self) -> np.ndarray:
        """

        :return: additional parameters like time etc.
        """
        return self.__params

    def set_params(self, p: np.ndarray):
        """

        :param p: additional parameters like time etc.
        :return: none
        """
        self.__params = p.copy()
