import numpy

from model.abstract_model import AbstractModel


class FreeFallModel(AbstractModel):

    def __init__(self):
        super().__init__(3)

    def model(self, a: list, b: list):
        return a[0] * b[0] * b[0] + a[1] * b[0] + a[2]

    def b_matrix(self, a: list, o: list):
        tmp = numpy.array((3, len(o)))
