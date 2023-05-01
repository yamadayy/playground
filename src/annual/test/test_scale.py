import math

from annual.scale import Scale


def test_init():
    scale = Scale(0.0, 1.0, 0.0)
    x = 1.0
    assert math.isclose(x, scale.observe(x), rel_tol=1e-10)


def test_init2():
    scale = Scale(0.1, 1.0, 0.0)
    x = 1.0
    assert not math.isclose(x, scale.observe(x), rel_tol=1e-10)
    assert math.isclose(x, scale.true_value(scale.observe(x)), rel_tol=1e-10)
