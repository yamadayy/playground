from abc import ABC, abstractmethod
from model.no_constraint import NoConstraint
import numpy as np


class AbstractModel(ABC):
    """
    This class is abstract class of the model which is described by a set of
    condition equations, f_i({x}_i, {a}_i) = 0. Let {x} be a set of unknown
    true value of {x}_0 which is obtained by direct observation. The vector
    v_i = x_i - x_0i is the vector of the observing errors. Assume that {a} is
    a set of parameters which we would like to estimate. A set {x}_i is subset
    of {x} and a set {a}_i is a subset of {a}. In this routine we assume simply
    that {x}_i = {x_i}. The solution of the least square problem consists in
    finding sets {x} and {a} which minimize the quadratic form v^T sigma^(-1) v
    while at the same time satisfying the condition equation, where sigma is
    covariance matrix of the observation.
    In this routine, the matrix X in the below paper which is defined by partial
    derivative of G with respect to x becomes unity from the assumtioon of
    {x}_i = {x_i}.
    See "Least Square Adjustment with relatively large observation errors,
    inaccurate initial approximations, or both" by Heinrich Eichhorn and Warren
    G. Clay, Mon. Not. R. astr. Soc. (1974) 166, 425-432.
    """
    def __init__(self, n: int):
        self.__n = n
        self.__constraint = NoConstraint()

    @abstractmethod
    def model(self, a0: np.ndarray, b: np.ndarray) -> np.ndarray:
        """
        The model of the above paper is describes as a number of condition
        equations, \vec{F}(\vec{a}, \vec{x}) = 0. This vector equation is
        separated into two parts, \vec{G}(\vec{a}, \vec{x}) = 0 and
        \vec{H}(\vec{a}) = 0. The condition equations which do not contain
        \vec{x} are denoted by \vec{H}, and other equations are denoted by
        \vec{G}. If we input \vec{a}_0 and \vec{x}_0 into the model,  \vec{G}_0
        and \vec{H}_0 are not zeros. This function returns the value of
        \vec{G}_0.
        :param a0: trial values of parameters which we would like to estimate
        :param b: additional parameters like time etc.
        :return: value(s) of the model, G_0.
        """
        pass

    @abstractmethod
    def b_matrix(self, a0: np.ndarray, o: np.ndarray) -> np.ndarray:
        """
        Matrix B is defined in the above paper as the partial derivative of G
        with respect to a.

        :param a0: trial values of parameters which we would like to estimate
        :param o: array of Observation objects.
        :return: Matrix B.
        """
        pass

    def num_param(self):
        """

        :return: number of parameters which we would like to estimate.
        """
        return self.__n

    def num_constraint(self):
        """
        This routine simply call the method with same name in ConcreteConstraint
         class.
        :return: number of constraints.
        """
        return self.__constraint.dimension()

    def c_matrix(self, a0: np.ndarray):
        """
        This routine simply call the method with same name in ConcreteConstraint
         class.
        :param a0: trial value of parameters which we would like to estimate.
        :return: constraint matrix C
        """
        return self.__constraint.c_matrix(a0)

    def constraint(self, a0: np.ndarray) -> np.ndarray:
        """
        This routine simply call the method with same name in ConcreteConstraint
         class.
        :param a0: trial value of parameters which we would like to estimate.
        :return: value(s) of constraints.  H_0
        """
        return self.__constraint.constraint(a0)
