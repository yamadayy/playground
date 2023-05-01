import math
import random
from annual.client import make_stars, make_scales, make_reference_stars


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
