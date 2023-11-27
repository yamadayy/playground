import numpy as np
from abc import ABC, abstractmethod

from jasmine_toolkit.optimize.abstract_model import AbstractModel
from jasmine_toolkit.optimize.observation import Observation


class AbstractSimulation(ABC):
    def __init__(self, model: AbstractModel, data_model):
        self.__model = model
        self.__data_model = data_model
        self.__obs = np.empty(0)

    def simulation(self, paramlist, error):
        n = paramlist.size
        for i in range(n):
            rr = np.random.normal(loc=0, scale=error)
            o = Observation(self.__model.model(self.__data_model, paramlist[i]) + rr, error, paramlist[i])
            self.__obs = np.append(self.__obs, np.array([o]))

    def get_observation_list(self):
        return self.__obs

    @abstractmethod
    def get_data_model(self):
        pass
