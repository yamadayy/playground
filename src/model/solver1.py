from model.abstract_model import AbstractModel
from model.abstract_solver import AbstractSolver


class Solver(AbstractSolver):
    def __init__(self, m: AbstractModel):
        super.__init__(m)