from jasmine_toolkit.converter.abstract_simulation import AbstractSimulation
from jasmine_toolkit.optimize.abstract_model import AbstractModel
from jasmine_toolkit.datamodel.toy.free_fall_model import FreeFallModel
import numpy as np


class GtoHeight(AbstractSimulation):
    def get_data_model(self):
        pass


if __name__ == '__main__':
    #  array a should be subclass of DataModel
    a = np.array([1.0, 0.0, 0.0])
    n = 11
    b = np.empty((n, 1))
    for i in range(n):
        b[i] = np.array([i * 0.1])
    error = 0.01

    m: AbstractModel = FreeFallModel()
    c = GtoHeight(m, a)
    c.simulation(b, error)
    oo = c.get_observation_list()

    for i in range(n):
        print(oo[i].get_value())
