class Scale:
    def __init__(self, o: float, s: float, t: float):
        self.__origin = o
        self.__scale = s
        self.__time = t

    def observe(self, pos: float) -> float:
        """
        Convert from true value to observed value
        :param pos: true position
        :return: observed value
        """
        return (pos - self.__origin) / self.__scale

    def true_value(self, observed: float) -> float:
        """
        Convert from observed value to true value.
        :param observed: observed value
        :return: true coordinate
        """
        return observed * self.__scale + self.__origin

    def get_time(self):
        return self.__time
