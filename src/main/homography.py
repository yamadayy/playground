import csv
import math
import numpy as np
import cv2


def affine(_p_original, _v_over_c, _v_y):
    p_trans = np.zeros(_p_original.shape, dtype="float32")
    for i in range(0, _p_original.shape[0]):
        p_trans[i] = relativistic_shift(_p_original[i], _v_over_c, _v_y)
    # print(p_trans)
    return cv2.getPerspectiveTransform(_p_original, p_trans)


def relativistic_shift_2(_original, _v_over_c, _vy):
    length = math.sqrt(_original[0] ** 2 + (_original[1] + _vy) ** 2)
    x_fac = _original[0] / length
    y_fac = (_original[1] + _vy) / length
    shift = _v_over_c * math.degrees(math.sin(math.radians(length)))
    p_trans = np.zeros(2, dtype="float32")
    p_trans[0] = _original[0] - x_fac * shift
    p_trans[1] = _original[1] - y_fac * shift
    return p_trans


def relativistic_shift(_original, _v_over_c, _vy):
    _x = math.radians(_original[0])
    _y = math.radians(_original[1])
    _theta_v = math.pi / 2 - math.radians(_vy)
    _theta = math.acos(math.sqrt(_x * _x + _y * _y))
    if _x == 0 and _y == 0:
        _phi = 0
    else:
        _phi = math.atan2(_y, -_x)
    _theta_c = math.atan2(math.sqrt(1 - _y * _y), math.sqrt(1 - _x * _x) * _y)
    _a = math.acos(math.cos(_theta_c) * math.cos(_theta) * math.sin(_phi) + math.sin(_theta_c) * math.sin(_theta))
    _b = _theta_c - _theta_v
    _psi = math.acos(math.cos(_theta_v) * math.cos(_theta) * math.sin(_phi) + math.sin(_theta_v) * math.sin(_theta))
    # _cos_B = (math.cos(_b) - math.cos(_psi) * math.cos(_a))/(math.sin(_psi) * math.sin(_a))
    if _a == 0:
        _cos_B = 0
    else:
        _cos_B = (math.cos(_psi) - math.cos(_a) * math.cos(_b)) / (math.sin(_a) * math.sin(_b))
    if _x == 0:
        _angle_B = 0
    elif _x > 0:
        _angle_B = math.acos(_cos_B)
    else:
        _angle_B = math.pi - math.acos(_cos_B)
    p_trans = np.zeros(2, dtype="float32")
    p_trans[0] = _original[0] - math.degrees(math.sin(_psi)) * _v_over_c * math.cos(_angle_B)
    p_trans[1] = _original[1] - math.degrees(math.sin(_psi)) * _v_over_c * math.sin(_angle_B)
    # p_trans[0] = math.degrees(math.sin(_psi)) * _v_over_c * math.cos(_angle_B)
    # p_trans[1] = math.degrees(math.sin(_psi)) * _v_over_c * math.sin(_angle_B)
    return p_trans


def homography(_p_original, _v_over_c, _v_y):
    A = np.zeros((2 * _p_original.shape[0], 8))
    b = np.zeros(2 * _p_original.shape[0])
    for i in range(0, _p_original.shape[0]):
        tmp = relativistic_shift(_p_original[i], _v_over_c, _v_y)
        b[2 * i] = tmp[0]
        b[2 * i + 1] = tmp[1]
    for i in range(0, _p_original.shape[0]):
        _x = _p_original[i, 0]
        _y = _p_original[i, 1]
        A[2 * i, 0] = 1
        A[2 * i, 1] = _x
        A[2 * i, 2] = _y
        A[2 * i, 3] = _x * _y
        A[2 * i + 1, 4] = 1
        A[2 * i + 1, 5] = _x
        A[2 * i + 1, 6] = _y
        A[2 * i + 1, 7] = _x * _y
    return np.linalg.pinv(A.T @ A) @ A.T @ b


def all_order(_p_original, _v_over_c, _v_y, _order):
    A = np.zeros((2 * _p_original.shape[0], (_order + 1) * (_order + 2)))
    b = np.zeros(2 * _p_original.shape[0])
    for i in range(0, _p_original.shape[0]):
        tmp = relativistic_shift(_p_original[i], _v_over_c, _v_y)
        b[2 * i] = tmp[0]
        b[2 * i + 1] = tmp[1]
    _s = int((_order + 1) * (_order + 2) / 2)
    for i in range(0, _p_original.shape[0]):
        _x = _p_original[i, 0]
        _y = _p_original[i, 1]
        A[2 * i, 0] = 1
        A[2 * i + 1, _s] = 1
        for j in range(1, _order + 1):
            for k in range(0, j + 1):
                A[2 * i, 2 * j - 1 + k] = _x ** (j - k) * _y ** k
                A[2 * i + 1, 2 * j - 1 + k + _s] = _x ** (j - k) * _y ** k
    # print(np.linalg.det(A.T @ A))
    return np.linalg.pinv(A.T @ A) @ A.T @ b


# M = affine(np.float32([[-0.25, 0.25], [0.25, 0.25], [0.25, -0.25], [-0.25, -0.25]]), 0.05, 4)
# print(M)
# p = homography(np.float32([[-0.25, 0.25], [0.25, 0.25], [0.25, -0.25], [-0.25, -0.25]]), 1e-4, 4)
# p = all_order(np.float32([[-0.25, 0.25], [0.25, 0.25], [0.25, -0.25], [-0.25, -0.25], [0, 0.25], [0, -0.25], [0, 0],
#                            [-0.2, -0.1], [-0.1, 0.2], [0.1, -0.], [0, 0.1]]), 0.05, 4, 2)
order = 4
trial = 100
n_component = (order + 1) * (order + 2)
n_points = n_component * 2
p = np.zeros((trial, n_component))
angle = 10
v_over_c = 1e-4
for i in range(0, trial):
    v = np.random.random((n_points, 2))
    vv = 0.5 * (v - 0.5)
    p[i] = all_order(vv, v_over_c, angle, order)
p_ave = np.mean(p, axis=0)
p_dev = np.std(p, axis=0)
# print(p_ave)
for i in range(0, n_component):
    if abs(p_ave[i]) < p_dev[i]:
        p_ave[i] = 0
print(p_ave)
print(p_dev)
# with open('a.csv', 'w') as f:
#     writer = csv.writer(f, delimiter=',')
#     writer.writerow(p)
# print(np.mean(a, axis=0))
# print(np.std(a, axis=0))
