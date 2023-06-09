from model.abstract_constraint import AbstractConstraint

class NoConstraint(AbstractConstraint):
    def constraint(self, a: list):
        return None

    def dimension(self):
        return 0

    def c_matrix(self, a: list):
        return None
