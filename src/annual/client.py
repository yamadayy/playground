import math
import random

import numpy as np

from annual.star import Star
from annual.scale import Scale


def estimate_scale_parameters(scale0: Scale, ref: list):
    sum_x = 0
    sum_xx = 0
    sum_y = 0
    sum_xy = 0
    n = 0
    for i in range(len(ref)):
        x = scale0.observe(ref[i].position()) + random.normalvariate(
            mu=0, sigma=observation_accuracy)
        y = ref[i].catalogue_position()
        sum_x = sum_x + x
        sum_xx = sum_xx + x * x
        sum_y = sum_y + y
        sum_xy = sum_xy + x * y
        n = n + 1
    mat = np.array([[sum_xx, sum_x], [sum_x, n]])
    vec = np.array([[sum_xy], [sum_y]])
    ans = np.linalg.inv(mat).dot(vec)
    return [ans[0][0], ans[1][0]]


def make_scales(n_season: int, n_exposure: int, dev_origin: float,
                dev_scale: float, dev_exposure: float):
    scales = []
    for i in range(n_season):
        t = i * 0.5 - (n_season - 1) * 0.25
        for j in range(n_exposure):
            dt = -0.125 + j * 0.25 / (n_exposure - 1.0)
            tt = t + dt
            origin = random.normalvariate(mu=0, sigma=dev_origin)
            scale = 1 + dev_scale * math.sin(2 * math.pi * tt)\
                    + random.normalvariate(mu=0, sigma=dev_exposure)
            scales.append(Scale(origin, scale, tt))
    return scales


def make_reference_stars(num_ref_stars: int, scale: float, ref_accuracy: float):
    refs = []
    for i in range(num_ref_stars):
        x = random.random() * scale
        y = random.normalvariate(mu=0, sigma=ref_accuracy)
        refs.append(Star(x, data_error=y))
    return refs


def make_stars(n: int, scale: float):
    targets = []
    for i in range(n):
        x = random.random() * scale
        targets.append(Star(x))
    return targets


def estimate_parallax(param: list, scales, s: Star):
    a = []
    for i in range(len(param)):
        a.append([1, param[i][2], math.sin(2 * math.pi * param[i][2])])
    mat_a = np.array(a)
    b = [[0],[0],[0]]
    for i in range(len(param)):
        observed_value = scales[i].observe(s.position())
        estimated_position = param[i][0] * observed_value + param[i][1]
        b[0][0] = b[0][0] + estimated_position
        b[1][0] = b[1][0] + estimated_position * param[i][2]
        b[2][0] = b[2][0] + estimated_position \
                  * math.sin(2 * math.pi * param[i][2])
    mat_b = np.array(b)
    astrometric_param = np.linalg.inv(mat_a.transpose().dot(mat_a)).dot(mat_b)
    return astrometric_param[2][0]


if __name__ == '__main__':
    num_star = 50
    num_reference_star = 10
    line_scale = 100
    amplitude_annual_scale = 0.25
    amplitude_every_observation = 0.01
    sigma_origin = 2
    num_exposure_per_season = 50
    num_season = 6
    observation_accuracy = 1
    reference_accuracy = 0.1
    target_accuracy = 0.03
    if False:
        observation_accuracy = 0
        reference_accuracy = 0
        amplitude_every_observation = 0
    mean = []
    stdev = []
    for rep in range(100):
        parallax = []
        stars = make_stars(num_star, line_scale)
        ref_stars = make_reference_stars(num_reference_star, line_scale,
                                         reference_accuracy)
        scales0 = make_scales(num_season, num_exposure_per_season, sigma_origin,
                              amplitude_annual_scale / line_scale,
                              amplitude_every_observation / line_scale)
        estimated_plate_param = []
        for k in range(len(scales0)):
            ans = estimate_scale_parameters(scales0[k], ref_stars)
            estimated_plate_param.append([ans[0], ans[1], scales0[k].get_time()])
        for i in range(len(stars)):
            parallax.append(estimate_parallax(estimated_plate_param, scales0,
                                              stars[i]))
        p = np.array(parallax)
        mean.append(p.mean())
        stdev.append(p.std())
    m = np.array(mean)
    s = np.array(stdev)
    print(str(m.mean()) + "," + str(s.mean()))
