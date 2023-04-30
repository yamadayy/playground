import random

import numpy as np

from star import Star
from scale import Scale

if __name__ == '__main__':
    num_star = 10
    num_reference_star = 10
    line_scale = 100
    amplitude_annual_scale = 0.25
    amplitude_every_observation = 0.01
    num_exposure_per_season = 50
    observation_accuracy = 1
    reference_accuracy = 0.1
    target_accuracy = 0.03
    stars = []
    for i in range(num_star):
        x = random.random() * line_scale
        # print(x)
        stars.append(Star(x))
    ref_stars = []
    for i in range(num_reference_star):
        x = random.random() * line_scale
        y = random.normalvariate(mu=0, sigma=reference_accuracy)
        # print(y)
        ref_stars.append(Star(x, data_error=y))

    sample = np.empty(0)
    for i in range(10000):
        sample = np.append(sample, random.normalvariate(mu=-1, sigma=0.3))
    print(sample.mean())
    print(sample.std())
