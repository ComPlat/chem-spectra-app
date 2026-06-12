import numpy as np
from chem_spectra.lib.shared.calc import get_curve_endpoint


def test_get_curve_endpoint_includes_trailing_points():
    xs = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    ys = np.array([0.0, 1.0, 2.0, 1.0, 0.0])
    iL, iU = get_curve_endpoint(xs, ys, 2.0, 4.0)
    assert iL == 2
    assert iU == 5
