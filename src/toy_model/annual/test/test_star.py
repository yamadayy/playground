import math

from toy_model.annual.star import Star


def test_star():
    x = 10
    s = Star(x)
    assert math.isclose(x, s.position(), rel_tol=1e-10)
    assert math.isclose(x, s.catalogue_position(), rel_tol=1e-10)


def test_star2():
    x = 10
    error = 0.1
    s = Star(x, data_error=error)
    assert math.isclose(x, s.position(), rel_tol=1e-10)
    assert math.isclose(x + error, s.catalogue_position(), rel_tol=1e-10)
