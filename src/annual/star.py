class Star:
    def __init__(self, x: float, data_error: float = 0):
        self.__position = x
        self.__data_error = data_error

    def position(self) -> float:
        return self.__position

    def catalogue_position(self) -> float:
        return self.__data_error + self.__position
