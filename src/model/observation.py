import pylab as p


class Observation:
    def __init__(self, v: float, s: float, p: list):
        self.__value: float = v
        self.__sigma: float = s
        self.__params = None
        if p != None:
            self.__params = p.copy()

    def get_value(self):
        return self.__value

    def get_sigma(self):
        return self.__sigma

    def get_params(self):
        return self.__params
