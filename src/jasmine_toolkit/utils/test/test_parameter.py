from jasmine_toolkit.utils import Parameters


def test_refer():
    p = Parameters()
    p.ready()
    assert p.n_ch == 16
