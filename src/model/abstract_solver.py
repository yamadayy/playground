from abc import ABC, abstractmethod
from model.abstract_model import AbstractModel
import numpy as np


class AbstractSolver(ABC):
    """
    This routine calculate the correction Delta of the trial set of {a}_0.
    See "Least Square Adjustment with relatively large observation errors,
    inaccurate initial approximations, or both" by Heinrich Eichhorn and Warren
    G. Clay, Mon. Not. R. astr. Soc. (1974) 166, 425-432.
    """
    def __init__(self, m: AbstractModel):
        self._model: AbstractModel = m
        self._np = self._model.num_param()
        self._nc = self._model.num_constraint()

    @abstractmethod
    def correction(self, a0: np.ndarray, o: np.ndarray):
        """
        This is the abstract method for returning correction.
        :param a0: trial values of parameters which we would like to estimate
        :param o: a set of observations
        :return: correction vector Delta.
        """
        pass

    @abstractmethod
    def covariant(self, a0: np.ndarray, o: np.ndarray):
        """
        This routine returns the covariant matrix of estiamted parameters.
        :param a0: trial values of parameters which we would like to estimate
        :param o: a set of observations
        :return: covariant matrix Sigma.
        """
        pass
