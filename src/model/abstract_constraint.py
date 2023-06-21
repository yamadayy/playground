from abc import ABC, abstractmethod
import numpy as np


class AbstractConstraint(ABC):
    """
    This is the abstraction of constraint.
    See "Least Square Adjustment with relatively large observation errors,
    inaccurate initial approximations, or both" by Heinrich Eichhorn and Warren
    G. Clay, Mon. Not. R. astr. Soc. (1974) 166, 425-432.
    """
    @abstractmethod
    def constraint(self, a0: np.ndarray) -> np.ndarray:
        """
        The constraints in the above paper are denoted by \vec{H}(\vec{a}) = 0.
        When we input to \vec{a}_0 (approximate value), H_0 = H(a_0) is not
        zero. This routine returns H_0 = H(a_0).
        :param a0: trial values of parameters which we would like to estimate
        :return: the values of constraint(s), H_0
        """
        pass

    @abstractmethod
    def dimension(self) -> int:
        """

        :return: number of constraints.
        """
        pass

    @abstractmethod
    def c_matrix(self, a0: np.ndarray) -> np.ndarray:
        """
        The matrix C in the above paper is defined as the partial derivative
        of H with respect to a.

        :param a0: trial values of parameters which we would like to estimate
        :return: Matrix C.
        """
        pass
