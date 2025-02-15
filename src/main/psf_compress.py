import math
import random
import numpy as np
from matplotlib import pyplot as plt


class StellarImage:
    def __init__(self, _nx, _ny):
        self._nx = _nx
        self._ny = _ny
        self.image = None

    def encode(self, image, base):
        self.image = image
        return 1

    def decode(self, stream, base):
        self.image = np.zeros((2 * _nx + 1) * (2 * _ny + 1))
        return 0


class StellarImageCollection:
    def __init__(self):
        self._data = None

    def generate_data(self, _n_sample):
        return None

    def pca(self):
        None


class GaussImageCollection(StellarImageCollection):
    def __init__(self):
        None


class ImageDataCollection(StellarImageCollection):
    def __init__(self):
        None


class Normal2d:
    def __init__(self, _stdev, _nx, _ny, _aspect=1):
        self.nx = _nx
        self.ny = _ny
        self.x0 = 0
        self.y0 = 0
        self.stdev = _stdev
        self.aspect = _aspect
        self.c0 = 1 / math.pi / 2 / _stdev / _stdev
        self.c1 = 1 / 2 / _stdev / _stdev
        self.X = None
        self.Y = None

    def set_meshgrid(self):
        _x = np.arange(-self.nx, self.nx + 0.1)
        _y = np.arange(-self.ny, self.ny + 0.1)
        self.X, self.Y = np.meshgrid(_x, _y, indexing='ij')

    def set_center_and_stdev(self, _x0, _y0, _c):
        self.x0 = _x0
        self.y0 = _y0
        self.stdev = _c

    def value(self, _x, _y):
        return self.c0 * math.exp(-((_x - self.x0) ** 2 + (_y - self.y0) ** 2) * self.c1)

    def get(self):
        _r = random.random()
        _t = random.random()
        x = math.sqrt(-2 * math.log(_r)) * self.stdev * math.cos(2 * math.pi * _t)
        y = math.sqrt(-2 * math.log(_r)) * self.stdev * math.sin(2 * math.pi * _t) * self.aspect
        return x + self.x0, y + self.y0

    def get_image(self, _n):
        _Z = np.zeros((2 * self.nx + 1, 2 * self.ny + 1))
        # read out 15e = 225 Poisson, dark 25x12.5=312, stray 5x12.5=62.5
        # expectation value is 600, poisson, but for read out, minus 210 (225 - 15)
        for i in range(2 * self.nx + 1):
            for j in range(2 * self.ny + 1):
                _Z[i][j] = 0  # np.random.poisson(600) - 210
        for i in range(_n):
            px, py = self.get()
            ix = round(px + self.nx)
            iy = round(py + self.ny)
            if (0 <= ix <= 2 * self.nx) and (0 <= iy <= 2 * self.ny):
                _Z[ix][iy] = _Z[ix][iy] + 1
        return _Z

    def get_meshgrid(self):
        return self.X, self.Y

    def get_flatten_size(self):
        return (2 * self.nx + 1) * (2 * self.ny + 1)


def pca(X):
    # X = X - X.mean(axis=0)
    cov = np.cov(X, rowvar=False)
    l, v = np.linalg.eig(cov)
    l_index = np.argsort(l)[::-1]
    v_ = v[:, l_index]
    return v_.T


def generate_data(_n, _normal):
    _normal.set_meshgrid()
    _data = np.empty((0, _normal.get_flatten_size()))
    for i in range(_n):
        _normal.set_center_and_stdev(random.random() - 0.5, random.random() - 0.5,
                                     (random.random() - 0.5) * 0.1 + 0.475)
        # maximum is 80000 electron = peak pixel contains 270000 photons,
        # but it should be 16000 (14bit -> total is 54000) or 65000 (16bit -> total is 220000)
        _n_photon = int(54000*0.01*pow(10, -0.4*math.log(random.random()+1/32, 2)))
        Z = _normal.get_image(_n_photon)
        _data = np.append(_data, Z.flatten().reshape(1, Z.flatten().shape[0]), axis=0)
    return _data


def show_contour(_x, _y, _z, _n):
    fig = plt.figure()
    for i in range(0, _n):
        ax = fig.add_subplot(2, int(_n / 2 + 0.5), i + 1)
        ax.set(aspect=1)
        ax.contour(_x, _y, _z[i].reshape((_x.shape[0], _x.shape[1])))
    plt.show()


def show_3d_graph(_x, _y, _z, _n):
    fig = plt.figure()
    for i in range(0, _n):
        ax = fig.add_subplot(2, int(_n / 2 + 0.5), i + 1, projection='3d')
        ax.plot_surface(_x, _y, _z[i].reshape((_x.shape[0], _x.shape[1])), linewidth=1)
    plt.show()


def rice_encode(_x, _k):
    m = 2 ** _k
    _r = _x % m  # if golomn, _x - q * m
    _q = _x >> _k
    _rem = bin(_r)[2:]
    _rem = '0' * (_k - len(_rem)) + _rem
    _quo = '0' * _q + '1'
    return _quo + _rem


def rice_decode(_str, _k):
    i = 0
    while _str[i] == '0':
        i = i + 1
    _r = 0
    for j in range(_k):
        _r = _r * 2
        if _str[j + i + 1] == '1':
            _r = _r + 1
    return _r + i * 2 ** _k


def calc_k(_x, _m):
    k = 0
    while (_m << k) < _x:
        k = k + 1
    return k


def show_principle_component(_nx, _ny):
    normal = Normal2d(0.6, _nx, _ny, _aspect=_ny/_nx)
    normal.set_meshgrid()
    XX, YY = normal.get_meshgrid()
    data = generate_data(_nx * _ny * 5, normal)
    T = pca(data)
    show_contour(XX, YY, T, 5)
    show_3d_graph(XX, YY, T, 5)


if __name__ == '__main__':
    _n_sample = 200
    _nx = 4
    _ny = 4
    _npix = (2 * _nx + 1) * (2 * _ny + 1)
    data = [1,2,3,4,5,6,7,8]
    x = np.abs(data).sum()
    n_data = len(data)
    k = calc_k(x, n_data)

"""
    normal = Normal2d(0.6, _nx, _ny, _aspect=1)
    data = generate_data(_n_sample, normal)  # shape = (_n_sample, _n_pix), 100 sample, 54000 photon->0.7 sec
    show_principle_component(5, 5)

    # print(data[1])

    mode = 4
    if mode == 1:
        # 最適なkの違い
        print("one image: npix={}, k={}".format(_npix, calc_k(np.abs(data[0]).sum(), _npix)))
        x = 0
        for i in range(_n_sample):
            x = x + np.abs(data[i]).sum()
        print("all image: npix={}, k={}".format(_npix, calc_k(x, _npix * _n_sample)))
    elif mode == 2:
        # kを変えるとどうなるか？
        for i in range(14):
            length = 0
            for j in range(45):
                length = length + len(rice_encode(int(data[1][j]), i))
            print("k = {}, {} bit".format(i, length))
    elif mode == 3:
        # 圧縮率のTの項数依存性
        T = pca(data)  # shape = (_n_pix, _n_pix)
        for n in range(45):
            kk_min = 100
            kk_max = 0
            coeff = np.zeros(n)
            res = np.zeros(45)
            length = 0
            for k in range(_n_sample):
                for i in range(n):
                    coeff[i] = np.dot(T[i], data[k])
                res = data[k]
                for i in range(n):
                    res = res - coeff[i] * T[i]
                _x = coeff.sum() + res.sum()
                kk = calc_k(_x, n + 45)
                if kk_min > kk:
                    kk_min = kk
                if kk_max < kk:
                    kk_max = kk
                for j in range(n):
                    length = length + len(rice_encode(int(coeff[j]), kk))
                for j in range(45):
                    length = length + len(rice_encode(int(res[j]), kk))
            print("n kmin kmax len  = ,{},{}, {}, {}".format(n, kk_min, kk_max, length / _n_sample))
    # 最適なkを使うように変更
    # nはTの展開次数、kはデータの番号、iはc[i]T[i]のi、jはn + 45の数字

    # show_principle_component(2, 4)
"""