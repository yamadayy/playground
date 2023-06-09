import math
import random
from annual.client import make_stars, make_scales, make_reference_stars,\
    estimate_parallax
from annual.scale import Scale
from annual.star import Star


def test_stars(monkeypatch):
    monkeypatch.setattr(random, 'random', lambda: 0.5)
    stars = make_stars(1, 10.0)
    assert len(stars) == 1
    assert math.isclose(stars[0].position(), 5.0, rel_tol=1e-10)
    assert math.isclose(stars[0].catalogue_position(), 5.0, rel_tol=1e-10)


def test_reference_stars(monkeypatch):
    monkeypatch.setattr(random, 'random', lambda: 0.5)
    monkeypatch.setattr(random, 'normalvariate', lambda mu=0, sigma=1: 0.1)
    stars = make_reference_stars(1, 10.0, 0.1)
    assert len(stars) == 1
    assert math.isclose(stars[0].position(), 5.0, rel_tol=1e-10)
    assert math.isclose(stars[0].catalogue_position(), 5.1, rel_tol=1e-10)


def test_make_scale(monkeypatch):
    monkeypatch.setattr(random, 'normalvariate', lambda mu=0, sigma=1: 0.1)
    scales = make_scales(1, 2, 1.0, 0.01, 1.0)
    assert len(scales) == 2
    assert math.isclose(scales[0].get_time(), -0.125, rel_tol=1e-7)
    assert math.isclose(scales[1].get_time(), 0.125, rel_tol=1e-7)
    assert math.isclose(scales[0].observe(5), 4.4833656, rel_tol=1e-7)
    assert math.isclose(scales[1].observe(5), 4.4260934, rel_tol=1e-7)


def test_estimate_parallax():
    estimated_scales = []
    for i in range(11):
        t = -0.5 + i * 0.1
        estimated_scales.append([1.0, 0.0, t])
    scales = []
    for i in range(len(estimated_scales)):
        scales.append(Scale(estimated_scales[i][1], estimated_scales[i][0],
                            estimated_scales[i][2]))
    star = Star(1.0)
    assert math.isclose(estimate_parallax(estimated_scales, scales, star), 0.0,
                        abs_tol=1e-7)