import numpy as np
from main.psf_compress import StellarImage, calc_k, encode_word, entropy


def test_calc_k():
    assert 1 == 1


def test_word_sequence():
    base = np.array([[0.42, 0.8, 0.42]])
    image = np.array([200, 400, 200])
    si = StellarImage(1, 1, base, image)
    ws = si.word_sequence()
    assert all(ws == np.array([488, -4, 9, -4]))


def test_calc_k():
    ws = np.array([488, -4, 9, -4])
    k = calc_k(ws)
    assert k == 7


def test_encode_word():
    assert encode_word(488, 7) == "011010001110"
    assert encode_word(424, 7) == "001010001110"
    assert encode_word(-9, 7) == "1000100110"


def test_encode():
    base = np.array([[0.42, 0.8, 0.42]])
    image = np.array([200, 400, 200])
    si = StellarImage(1, 1, base, image)
    assert si.encode() == "0111011010001110100001001000001001101000010010"


def test_entropy():
    assert entropy(np.array([1, 2, 3, 3])) == 1.5
