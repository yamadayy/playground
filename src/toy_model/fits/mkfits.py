import random

import numpy as np

nx = 100
ny = 100
n_stars = 10
max_photons = 10000
psf_width = 1.0
array = np.zeros((nx, ny), dtype='float64')

for s in range(n_stars):
    x = random.random() * nx
    y = random.random() * ny
    n_photons = int(random.random() * max_photons)
    for p in range(n_photons):
        px = random.gauss(x, psf_width)
        py = random.gauss(y, psf_width)
        if (0 < px < nx) and (0 < py < ny):
            array[int(px), int(py)] += 1

np.savetxt("tmp.csv", array, delimiter=',')
