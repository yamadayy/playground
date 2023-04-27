import random

from star import Star
from scale import Scale

if __name__ == '__main__':
    num_star = 10
    num_reference_star = 10
    line_scale = 100
    amplitude_annual_scale = 0.08
    amplitude_every_observation = 0.01
    num_exposure = 10
    target_accuracy = 0.03
    stars = []
    for i in range(num_star):
        x = random.random() * line_scale
        print(x)
        stars.append(Star(x))

    print()
    for i in range(len(stars)):
        print(stars[i].position())