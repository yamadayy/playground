import numpy as np
import math


def root_n_law(n_repeat: int, sequence: list):
    a = np.random.randn(n_repeat)
    print(str(a.mean()) + ", " + str(a.std()))
    for i in sequence:
        a = np.empty(0)
        for j in range(n_repeat):
            a = np.append(a, np.random.randn(i).mean())
        print(str(a.mean()) + ",  " + str(a.std()))


def model(n_repeat: int, sequence: list, v0: float):
    for i in sequence:
        x = np.empty(0)
        for j in range(n_repeat):
            t = np.linspace(0, 2, i)
            r = np.random.randn(i)
            h = 0.5 * t * t + v0 * t + 0.01 * r
            a = np.sum(h * t * t) / np.sum(t * t * t * t)
            x = np.append(x, a * t[int(i / 2)] * t[int(i / 2)])
        print(str(i) + ", " + str(math.fabs(x.mean() - 0.5)) + ", " + str(x.std()))


def model2(n_repeat: int, sequence: list, v0: float):
    for i in sequence:
        x = np.empty(0)
        hh = np.empty(0)
        for j in range(n_repeat):
            t = np.linspace(0, 2, i)
            r = np.random.randn(i)
            h = 0.5 * t * t + v0 * t + 0.01 * r
            a = np.sum(h * t * t) / np.sum(t * t * t * t)
            h0 = a * t * t - h
            x = np.append(x, a * t[int(i / 2)] * t[int(i / 2)])
            hh = np.append(hh, h0.mean())
        print(str(i) + ", " + str(math.fabs(hh.mean())) + ", " + str(x.std()))


if __name__ == '__main__':
    # root_n_law(100, [10, 30, 100, 300, 1000, 3000, 10000])
    v_ini = 0
    print("v0=" + str(v_ini) + ",compare with true value")
    model(100, [11, 31, 101, 301, 1001, 3001, 10001], v_ini)
    print("          compare with observed value")
    model2(100, [11, 31, 101, 301, 1001, 3001, 10001], v_ini)
    v_ini = 0.001
    print("v0=" + str(v_ini) + ",compare with true value")
    model(100, [11, 31, 101, 301, 1001, 3001, 10001], v_ini)
    print("          compare with observed value")
    model2(100, [11, 31, 101, 301, 1001, 3001, 10001], v_ini)
