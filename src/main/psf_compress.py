import math
import random
import matplotlib.pyplot as plt
import numpy as np


def calc_k(_data):
    """
    calculate k for Golomb Rice coding
    :param _data: numpy.ndarray of input data
    :return: the value of k
    """
    _a = np.abs(_data).sum()
    _n_data = len(_data)
    _k = 0
    while (_n_data << _k) < _a:
        _k = _k + 1
    return _k


def encode_word(_word, _k):
    """
    Golomb rice encoding from value to bit stream. Stream contains only '1' and '0' byte stream and should be
     reinterpreted as bit stream.
    :param _word: data
    :param _k: value of k
    :return: output stream
    """
    SYMBOL = '1'
    STOP = '0'
    format_string = "0{}b".format(_k)
    if _word > 0:
        _byte = format(_word, format_string)
        _lower = "0" + _byte[-_k:]
    else:
        _byte = format(-_word, format_string)
        _lower = "1" + _byte[-_k:]
    _upper = _byte[:-_k]
    _unary = ""
    for i in range(len(_upper) + 1):
        _unary += SYMBOL
    _unary += STOP
    return _lower + _unary


def encode_word_original(_word, _k):
    """
    The same as function encode_word which we assume all word is positive.
    :param _word:
    :param _k:
    :return:
    """
    SYMBOL = '1'
    STOP = '0'
    format_string = "0{}b".format(_k)
    _byte = format(_word, format_string)
    _lower = _byte[-_k:]
    _upper = _byte[:-_k]
    _unary = ""
    for i in range(len(_upper) + 1):
        _unary += SYMBOL
    _unary += STOP
    return _lower + _unary


class StellarImage:
    """
    Class which process one stellar image.
    """
    def __init__(self, _nx, _ny, _base, _image=None):
        """
        Constructor

        :param _nx: number of pixels in x direction
        :param _ny: number of pixels in y direction.
        :param _base: dimension 2 array. Eigen vectors of Karhunen-Loeve expansion basis.
        :param _image: serialized image values. Array size should be _nx * _ny
        """
        self._nx = _nx
        self._ny = _ny
        self._image = _image
        self._base = _base

    def word_sequence(self):
        """
        Convert pixel values of image to model parameters and residuals.

        :return: word_sequence for compression
        """
        c = np.empty(0)
        for i in range(self._base.shape[0]):
            c = np.append(c, np.dot(self._base[i], self._image))
        res = self._image
        for i in range(self._base.shape[0]):
            res = res - c[i] * self._base[i]
        c = np.append(c, res)
        c = c.astype(int)
        return c

    def encode(self):
        _ws = self.word_sequence()
        _k = calc_k(_ws)
        encoded_string = format(_k, "04b")
        for i in range(len(_ws)):
            encoded_string += encode_word(_ws[i], _k)
        return encoded_string

    def decode(self, stream):
        self._image = np.zeros((2 * self._nx + 1) * (2 * self._ny + 1))
        return 0


class StellarImageCollection:
    def __init__(self, _nx, _ny):
        """
        Super class of the stellar image correction.
        :param _nx: number of pixels in the stellar window in x direction
        :param _ny: number of pixels in the stellar window in y direction
        """
        self._data = None
        self.X = None
        self.Y = None
        self.nx = _nx
        self.ny = _ny

    def generate_data(self, _n_sample):
        return None

    def pca(self):
        """
        Principal Component Analysis.
        :return: Eigenvectors of the variance-covariance matrix
        """
        # X = X - X.mean(axis=0)
        cov = np.cov(self._data, rowvar=False)
        l, v = np.linalg.eig(cov)
        l_index = np.argsort(l)[::-1]
        v_ = v[:, l_index]
        return v_.T

    def show_pca(self, file_name):
        """
        Draw principal components
        :param file_name: csv file name for output
        :return: None
        """
        v = self.pca()
        vv = np.empty(0)
        for i in range(self.nx * self.ny):
            vv = np.append(vv, v[i])
        vv = vv.reshape(self.nx * self.ny, 7)
        with open(file_name, 'a') as f_handle:
            f_handle.truncate(0)
            for i in range(len(vv)):
                np.savetxt(f_handle, vv[i], delimiter=',')


class GaussImageCollection(StellarImageCollection):
    """
    class for make and manage gauss psf data
    """
    def __init__(self, _nx, _ny, _stddev, _aspect):
        """
        Constructor
        :param _nx: stellar window size in x direction
        :param _ny: stellar window size in y direction
        :param _stddev: standard deviation of gauss function in the unit of pixel number
        :param _aspect: aspect ratio of the stellar window
        """
        super().__init__(_nx, _ny)
        self.x0 = 0
        self.y0 = 0
        self.stddev = _stddev
        self.aspect = _aspect

    def _set_meshgrid(self):
        _x = np.arange(-int(self.nx / 2), int(self.nx / 2) + 0.1)
        _y = np.arange(-int(self.ny / 2), int(self.ny / 2) + 0.1)
        self.X, self.Y = np.meshgrid(_x, _y, indexing='ij')

    def _get(self):
        _r = random.random()
        _t = random.random()
        x = math.sqrt(-2 * math.log(_r)) * self.stddev * math.cos(2 * math.pi * _t)
        y = math.sqrt(-2 * math.log(_r)) * self.stddev * math.sin(2 * math.pi * _t) * self.aspect
        return x + self.x0, y + self.y0

    def _set_center_and_stddev(self, _x0, _y0, _c):
        self.x0 = _x0
        self.y0 = _y0
        self.stddev = _c

    def _get_image(self, _n):
        _Z = np.zeros((self.nx, self.ny))
        # read out 15e = 225 Poisson, dark 25x12.5=312, stray 5x12.5=62.5
        # expectation value is 600, poisson, but for read out, minus 210 (225 - 15)
        for i in range(self.nx ):
            for j in range(self.ny):
                _Z[i][j] = 0  # np.random.poisson(600) - 210
        for i in range(_n):
            px, py = self._get()
            ix = round(px + int(self.nx / 2))
            iy = round(py + int(self.ny / 2))
            if (0 <= ix <= self.nx - 1) and (0 <= iy <= self.ny - 1):
                _Z[ix][iy] = _Z[ix][iy] + 1
        return _Z

    def generate_data(self, _n_sample):
        """
        generate gauss psf window data
        :param _n_sample: number of generated samples
        :return: generated data of numpy.ndarray with shape (nx, ny, n_sample)
        """
        self._set_meshgrid()
        self._data = np.empty((0, self.nx * self.ny))
        # maximum is 80000 electron in peak pixel, image contains 270000 photons,
        # but it should be 16000 (14bit -> total is 54000) or 65000 (16bit -> total is 220000)
        _num_max_photon = 54000
        _mag_width = 7
        _coefficient_b = 1 / math.pow(2.0, _mag_width)
        _coefficient_a = math.pow(10, -0.4*_mag_width)
        for i in range(_n_sample):
            self._set_center_and_stddev(random.random() - 0.5, random.random() - 0.5,
                                        (random.random() - 0.5) * 0.03 + 0.475)
            _n_photon = int(_num_max_photon * _coefficient_a * pow(10, -0.4 * math.log(random.random()
                                                                                       + _coefficient_b, 2)))
            Z = self._get_image(_n_photon)
            self._data = np.append(self._data, Z.flatten().reshape(1, Z.flatten().shape[0]), axis=0)
        return self._data


class ImageDataCollection(StellarImageCollection):
    """
    generate data array from the csv file.
    """
    def __init__(self, _nx, _ny):
        """
        Constructor
        :param _nx: stellar window size in x direction
        :param _ny: stellar window size in y direction
        """
        super().__init__(_nx, _ny)

    def generate_data(self, _n_sample):
        """

        :param _n_sample:
        :return: numpy.ndarray with shapte (nn_sample, nx * ny)
        """
        a = np.loadtxt("data.csv", delimiter=',')
        line_count = len(a)
        self._data = a.flatten().reshape((int(line_count / 9), 99))
        return self._data

    def set_data(self, data):
        self._data = data


def show_statistics(data, v, half_nx):
    """
    show statistics of compression of various number of bases.
    :param data: data array
    :param v: Karhunen Loeve bases
    :param half_nx: half of pixel size in x direction
    :return: non
    """
    original_data_size = data.shape[0] * data.shape[1] * 16
    half_ny = int(data.shape[1] / half_nx)
    print(original_data_size)
    print("number of base, length")
    for i in range(len(v)):
        compress_length = 0
        for j in range(data.shape[0]):
            si = StellarImage(half_nx, half_ny, v[0:i], data[i])
            compress_length += len(si.encode())
        print("{}, {}".format(i, compress_length))


def gauss_psf_statistics(_nx, _ny, n_sample):
    gi = GaussImageCollection(_nx, _ny, 2.5, 2)
    data = gi.generate_data(n_sample)
    v = gi.pca()
    show_statistics(data, v, int(_nx / 2))


def data_statistics(_nx, _ny):
    idc = ImageDataCollection(_nx, _ny)
    data = idc.generate_data(-1)
    v = idc.pca()
    show_statistics(data, v, int(_nx / 2))


def selected_data_statistics():
    half_nx = 5
    half_ny = 4
    idc = ImageDataCollection(2 * half_nx + 1, 2 * half_ny + 1)
    data = idc.generate_data(-1)
    num_pix = (2 * half_nx + 1) * (2 * half_ny + 1)
    center_pixel_id = int(num_pix / 2)
    selected_data = np.empty(0)
    count_max = 0
    count_min = 10000000
    for i in range(len(data)):
        # print(data[i].argmax())
        if data[i].argmax() == center_pixel_id:
            total = np.sum(data[i])
            if count_max < total:
                count_max = total
            if count_min > total:
                count_min = total
            # data[i] = np.array(data[i], dtype=np.int32) >> 6
            selected_data = np.append(selected_data, data[i])
    selected_data = selected_data.flatten().reshape((int(len(selected_data) / 99), 99))
    idc.set_data(selected_data)
    v = idc.pca()
    print("max={}, min={}".format(count_max, count_min))
    show_statistics(selected_data, v, half_nx)


def statistics_with_common_k(data, step):
    """
    Statistics which we assume that the value of k is common for all data
    :param data: data
    :param step: selected steps
    :return: None
    """
    selected_data = data[::step].flatten()
    k = calc_k(selected_data)
    length = 0
    for i in range(len(selected_data)):
        encoded_string = encode_word_original(int(selected_data[i]), k)
        length += len(encoded_string)
    print("original data: {}, encoded bits: {}, average bits={}".format(i + 1, length, length / (i + 1)))


def entropy(_data):
    """
    Calculate and print information entropy from all data.
    :param _data: data numpy.ndarray
    :return: entropy
    """
    data = _data
    u, count = np.unique(data, return_counts=True)
    probability = count / data.size
    ent = np.sum(-probability * np.log2(probability))
    return ent


def show_probability_histgram(_data):
    """
    Show histgram of probability distribution
    :param _data: raw data
    :return: non
    """
    data = _data
    u, count = np.unique(data, return_counts=True)
    fig, ax = plt.subplots()
    ax.plot(u, np.log10(count))
    ax.set_title('distribution')
    plt.xlabel('data')
    plt.ylabel('log count')
    plt.show()


if __name__ == '__main__':
    # gauss_psf_statistics(9, 5, 1000)
    # data_statistics(11, 9)
    # selected_data_statistics()
    # print(entropy(np.loadtxt("data.csv", delimiter=',').flatten()))
    # show_probability_histgram(np.loadtxt("data.csv", delimiter=',').flatten())
    statistics_with_common_k(ImageDataCollection(11, 9).generate_data(-1), 1)
