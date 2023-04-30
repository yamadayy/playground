import math
import random

import numpy as np

from star import Star
from scale import Scale


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


def make_scales():
    scales = []
    for i in range(num_season):
        t = i * 0.5 - (num_season - 1) * 0.25
        for j in range(num_exposure_per_season):
            dt = -0.125 + j * 0.25 / (num_exposure_per_season - 1.0)
            tt = t + dt
            origin = random.normalvariate(mu=0, sigma=2)
            scale = 1 + amplitude_annual_scale * math.sin(2 * math.pi * tt) \
                    / line_scale
            scales.append(Scale(origin, scale))
    return scales


def make_reference_stars():
    refs = []
    for i in range(num_reference_star):
        x = random.random() * line_scale
        y = random.normalvariate(mu=0, sigma=reference_accuracy)
        # print(y)
        refs.append(Star(x, data_error=y))
    return refs


def make_stars():
    targets = []
    for i in range(num_star):
        x = random.random() * line_scale
        # print(x)
        targets.append(Star(x))
    return targets


if __name__ == '__main__':
    num_star = 10
    num_reference_star = 10
    line_scale = 100
    amplitude_annual_scale = 0.25
    amplitude_every_observation = 0.01
    num_exposure_per_season = 50
    num_season = 6
    observation_accuracy = 1
    reference_accuracy = 0.1
    target_accuracy = 0.03
    if True:
        observation_accuracy = 0
        reference_accuracy = 0
    stars = make_stars()
    ref_stars = make_reference_stars()
    scales = make_scales()
    print([scales[0]._Scale__scale, scales[0]._Scale__origin])
    ans = estimate_scale_parameters(scales[0], ref_stars)
    print(ans)
