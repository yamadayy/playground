import numpy as np


def root_n_law(n_repeat: int, sequence: list):
    a = np.random.randn(n_repeat)
    print(str(a.mean()) + ", " + str(a.std()))
    for i in sequence:
        a = np.empty(0)
        for j in range(n_repeat):
            a = np.append(a, np.random.randn(i).mean())
        print(str(a.mean()) + ", " + str(a.std()))


if __name__ == '__main__':
    root_n_law(100, [10, 30, 100, 300, 1000, 3000, 10000])