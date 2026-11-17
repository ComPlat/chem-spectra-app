"""Writing measured values as text, as short as their precision allows."""

import numpy as np


def _shortest(values):
    """Each value as the shortest decimal that keeps its precision.

    mzML and netCDF store m/z and intensity as float32, which reach us
    widened to float64. Printed as float64 they carry 17 significant digits
    (`100.18072509765625`) of which only the first 7 to 9 were ever measured.
    When every value is exactly a float32, write the shortest decimal that
    parses back to the same float32 (`100.180725`); read as float64 it lies
    within half a float32 step of the original. Anything else keeps Python's
    float64 repr, as before.

    numpy prints large float32 values in scientific notation
    (`1.2345679e+08`); those are rewritten in Python's repr style, which
    keeps positional notation below 1e16, so the text looks as it did.
    """
    values = np.asarray(values, dtype=np.float64)
    narrow = values.astype(np.float32)
    if not np.array_equal(narrow.astype(np.float64), values, equal_nan=True):
        return [repr(v) for v in values.tolist()]
    text = narrow.astype(str)
    scientific = np.char.find(text, 'e') >= 0
    if scientific.any():
        text = text.astype(object)
        text[scientific] = [repr(float(t)) for t in text[scientific]]
    return text.tolist()


def xy_lines(xs, ys):
    """`x, y` lines, one per point, for an `(XY..XY)` data table."""
    return [
        '{}, {}\n'.format(x, y)
        for x, y in zip(_shortest(xs), _shortest(ys))
    ]
