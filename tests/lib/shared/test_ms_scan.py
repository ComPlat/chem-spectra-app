"""scan_ratios gives exactly what the per-point loop it replaced gave."""
import numpy as np
import pytest

from chem_spectra.lib.shared.ms_scan import scan_ratios


def _loop(spc, low, high):
    # the loop as it was in the mzML, JCAMP and netCDF MS converters
    base, seed, beyond = [], [], []
    for pk in spc:
        if low < pk[0] < high:
            seed.append(pk[1])
        elif low + 1 < pk[0] < high + 1:
            seed.append(pk[1])
        elif low + 23 < pk[0] < high + 23:
            seed.append(pk[1])
        elif low + 39 < pk[0] < high + 39:
            seed.append(pk[1])
        if pk[0] <= high + 39:
            base.append(pk[1])
        elif high + 39 < pk[0]:
            beyond.append(pk[1])
    max_base = max(base, default=0.1)
    max_seed = max(seed, default=0)
    max_oorg = max(beyond, default=0)
    return 100 * max_seed / max_base, 100 * max_oorg / max_base, max_seed


@pytest.mark.parametrize('dtype', [np.float32, np.float64])
@pytest.mark.parametrize('mass', [150.0, 230.079907196, 301.1])
def test_same_numbers_and_types_as_the_loop(dtype, mass):
    rng = np.random.default_rng(1)
    for _ in range(50):
        n = int(rng.integers(0, 400))
        xs = rng.uniform(100, 400, n)
        ys = rng.uniform(0, 1e6, n)
        spc = np.column_stack([xs, ys]).astype(dtype)
        want = _loop(spc, mass - 0.5, mass + 0.5)
        got = scan_ratios(spc, mass - 0.5, mass + 0.5)
        assert got == want
        assert [type(v) for v in got] == [type(v) for v in want]


def test_an_empty_scan():
    assert scan_ratios(np.array([]), 229.5, 230.5) == (0.0, 0.0, 0)


def test_the_target_windows():
    # the ion, +1, +23 (Na), +39 (K) count as the target; +2 does not
    spc = np.array([[230.0, 10.0], [231.0, 20.0], [232.0, 99.0],
                    [253.0, 30.0], [269.0, 40.0], [300.0, 5.0]])
    ratio, noise_ratio, max_seed = scan_ratios(spc, 229.5, 230.5)
    assert max_seed == 40.0
    assert ratio == 100 * 40.0 / 99.0
    assert noise_ratio == 100 * 5.0 / 99.0
