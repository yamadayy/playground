import numpy as np

from model.abstract_model import AbstractModel


class FreeFallModel(AbstractModel):

    def __init__(self):
        super().__init__(3)

    def model(self, a: list, b: list):
        return a[0] * b[0] * b[0] + a[1] * b[0] + a[2]

    def b_matrix(self, a: list, o: list):
        tmp = np.zeros((len(o), 3))
        for i in range(len(o)):
            tmp[i][0] = o[i] * o[i]
            tmp[i][1] = o[i]
            tmp[i][2] = 1.
        return tmp


if __name__ == '__main__':
    model: AbstractModel = FreeFallModel()
    print('number of parameters:  ' + str(model.num_param()))
    print('number of constraints: ' + str(model.num_constraint()))
    n = 11
    oo = []
    aa = [1.0, 0.0, 0.0]
    error = 0.01
    r = np.random.normal(loc=0, scale=error, size=n)
    for ii in range(n):
        t = 0.1 * ii
        oo.append(model.model(aa, [t]) + 0 * r[ii])
    print('observed values are ' + str(oo))
    mat_b = model.b_matrix(aa, oo)
    print('B = ' + str(mat_b))
    print('C = ' + str(model.c_matrix(aa)))
    a_try = [0.8, 0.0, 0.0]
    sigma = np.identity(n) * error * error
    F0 = np.zeros((n, 1))
    for j in range(10):
        for ii in range(n):
            F0[ii][0] = model.model(a_try, [0.1 * ii]) - oo[ii]
        sigma_param = np.linalg.inv(np.dot(mat_b.T, np.dot(np.linalg.inv(sigma),
                                                           mat_b)))
        delta = - np.dot(sigma_param, np.dot(mat_b.T,
                                             np.dot(np.linalg.inv(sigma), F0)))
        for ii in range(3):
            a_try[ii] = a_try[ii] + delta[ii][0]
        print(str(delta) + " : " + str(a_try))
