import math
import random
import numpy as np


def func(xxxx: float) -> float:
    a = 1.0
    b = 0.1
    c: float = 0.1
    return a * xxxx * xxxx + b * xxxx + c


def func1(coeff: np.matrix, xx: float) -> float:
    return coeff[0][0] * xx * xx + coeff[1][0] * xx + coeff[2][0]


def generate_data(xxx: float) -> float:
    r = random.normalvariate(0, 0.01)
    return r + func(xxx)


def one_trial():
    p = 0
    q = 0
    r = 0
    d = 0
    e = 0
    f = 0
    g = 0
    h = 0
    for i in range(-10, 11):
        x1 = i * 0.1
        y = generate_data(x1)
        p = p + x1 * x1 * x1 * x1
        q = q + x1 * x1 * x1
        r = r + x1 * x1
        d = d + x1
        e = e + 1
        f = f + x1 * x1 * y
        g = g + x1 * y
        h = h + y
    mat_a = np.array([[p, q, r], [q, r, d], [r, d, e]])
    vec_b = np.array([[f], [g], [h]])
    vec_c = np.linalg.inv(mat_a).dot(vec_b)
    return vec_c


if __name__ == '__main__':
    sum1: float = 0
    sum2: float = 0
    n: float = 0
    x0 = 1
    true_value = func(x0)
    for j in [10, 30, 100, 300, 1000, 3000]:
        for i in range(j):
            coeff = one_trial()
            x = func1(coeff, x0)
            sum1 = sum1 + x
            sum2 = sum2 + x * x
            n = n + 1
        ave = sum1 / n
        stdev = math.sqrt(sum2 / n - ave * ave)
        print(str(true_value - ave) + "," + str(ave) + "," + str(stdev))
