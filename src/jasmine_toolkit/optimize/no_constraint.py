from jasmine_toolkit.optimize.abstract_constraint import AbstractConstraint
import numpy as np


class NoConstraint(AbstractConstraint):
    """ ConcreteConstraint class without constraint.

    This class is concrete subclass of AbstractConstraint where there is no
    constraint, i.e. the dimension of H_0 is zero.
    """
    def constraint(self, a0: np.ndarray) -> np.ndarray:
        """ constraint method for the class without constraint.

        :param a0: trial values of parameters which we would like to estimate
        :return: None
        """
        return None

    def dimension(self) -> int:
        return 0

    def c_matrix(self, a0: np.ndarray) -> np.ndarray:
        return None
