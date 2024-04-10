import math

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv


def airy(eps: float, _x: float):
    if abs(_x) < 1e-10:
        tmp = 1 - eps ** 2
    else:
        tmp = 2 * jv(1, _x) / _x - 2 * eps * jv(1, eps * _x) / _x
    return tmp * tmp


def integ_1d(x_min:float, x_max:float, eps:float, n:int):
    tmp = 0
    x = np.linspace(x_min, x_max, n)
    dx = (x_max - x_min) / n
    for xx in x:
        tmp = tmp + dx * airy(eps, xx)
    return tmp


def integ_r(x_max:float, eps:float, n:int):
    tmp = 0
    x_min = 0
    x = np.linspace(x_min, x_max, n)
    dx = (x_max - x_min) / n
    for xx in x:
        tmp = tmp + dx * airy(eps, xx) * 2 * math.pi * xx
    return tmp / (1 - eps ** 2) / 4 / math.pi


def integ_2d(x_min:float, x_max:float, y_min:float, y_max:float, eps:float, n:int):
    x = np.linspace(x_min, x_max, n)
    y = np.linspace(y_min, y_max, n)
    tmp = 0
    dv = (x_max - x_min) / n * (y_max - y_min) / n
    for xx in x:
        for yy in y:
            r = math.sqrt(xx * xx + yy * yy)
            tmp = tmp + dv * airy(eps, r)
    return tmp / 4 / math.pi / (1 - eps * eps)


def show_graph(x_min:float, x_max:float, eps:float):
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)
    ax.set_title("Bessel functions of the first kind", size=15)
    ax.grid()
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-1.5, 1.5)
    ax.set_xlabel("x", size=15, labelpad=10)
    ax.set_ylabel("Jv(x)", size=15, labelpad=8)
    x = np.linspace(x_min, x_max, 65)
    y = []
    for xx in x:
        y.append(jv(1,xx))
#       y.append(airy(eps, xx))
    ax.plot(x, y)
    plt.show()


# x_max = 0.9953
# print("x_max=" + str(x_max) + " integ=" + str(integ_2d(-x_max, x_max, -x_max, x_max, 0.35, 100)))
a = [[1,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15],[16,17,18,19,20],[21,22,23,24,25]]
print(a)
print(a[1:2])
print(np.sum(a[1:2]))