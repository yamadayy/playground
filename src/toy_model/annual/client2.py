import math
import random

import matplotlib.pyplot as plt
import numpy as np


def generate_stars(num: int):
    stars = []
    for i in range(num):
        stars.append(random.random() * 2 * ll - ll)
    return stars


def calculate_a_and_f(stars: list):
    size = len(stars)
    mat_a = np.empty((size, 3))
    vec_f = np.empty((size, 1))
    for i in range(size):
        true_value = a + b * stars[i] + c * stars[i] * stars[i]
        observed_value = random.normalvariate(mu=true_value, sigma=1)
        mat_a[i][0] = 1
        mat_a[i][1] = stars[i]
        mat_a[i][2] = stars[i] * stars[i]
        vec_f[i][0] = - observed_value
    return mat_a, vec_f


def observe(stars):
    size = len(stars)
    obs = []
    for i in range(size):
        true_value = a + b * stars[i] + c * stars[i] * stars[i]
        obs.append(random.normalvariate(mu=true_value, sigma=1))
    return obs


def est_error(a, observed, stars):
    true_value = []
    for i in range(len(observed)):
        true_value.append((math.sqrt(a[1] * a[1] - 4 * a[2] * (a[0] - observed[i]))
                           - a[1]) / 2 / a[2] - stars[i])
    return true_value


if __name__ == '__main__':
    ll = 100
    a = 3
    b = 1.01
    c = 1E-4
    m = 5
    k = 5
    n = 5
    gaia_stars = generate_stars(m)
    normal_stars = generate_stars(k)
    print("Gaia stars  :" + str(gaia_stars))
    print("Normal stars:" + str(normal_stars))
    stat = np.empty((k, n))
    for j in range(n):
        mat_a_gaia, vec_f_gaia = calculate_a_and_f(gaia_stars)
        cova_param = np.linalg.inv(np.transpose(mat_a_gaia).dot(mat_a_gaia))
        delta = cova_param.dot(np.transpose(mat_a_gaia).dot(vec_f_gaia))
        a_estimated = - delta[0][0]
        b_estimated = - delta[1][0]
        c_estimated = - delta[2][0]
        mat_a_normal = np.empty((m, 3))
        vec_f_normal = np.empty((m, 1))
        obs = observe(normal_stars)
        pos_error = est_error([a_estimated, b_estimated, c_estimated], obs, normal_stars)
        for i in range(k):
            stat[i][j] = pos_error[i]
    # print(stat)
    print(str(np.average(stat[0])) + ", " + str(np.std(stat[0])))
    plt.hist(stat[0])
    plt.show()
    print(str(np.average(stat[1])) + ", " + str(np.std(stat[1])))
    plt.hist(stat[1])
    plt.show()
    print(str(np.average(stat[2])) + ", " + str(np.std(stat[2])))
    plt.hist(stat[2])
    plt.show()
    print(str(np.average(stat[3])) + ", " + str(np.std(stat[3])))
    plt.hist(stat[3])
    plt.show()
    print(str(np.average(stat[4])) + ", " + str(np.std(stat[4])))
    plt.hist(stat[4])
    plt.show()
