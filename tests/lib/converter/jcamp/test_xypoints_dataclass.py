import numpy as np
import pytest

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter

# Deliberately non-uniform, non-monotonic X — a miniature cyclic sweep. An
# evenly-spaced reconstruction from FIRSTX/LASTX cannot reproduce it, which is
# what makes it a usable probe.
XS = [0.0, 0.1, 0.3, 0.6, 0.4, 0.2, 0.05]
YS = [1.0, 1.4, 2.1, 3.0, 2.4, 1.6, 1.1]


def _jcamp(data_class, table_label):
    rows = ''.join(
        '{:.6E}, {:.6E}\n'.format(x, y) for x, y in zip(XS, YS)
    )
    header = '' if data_class is None else '##DATA CLASS={}\n'.format(data_class)
    return (
        '##TITLE=Spectrum\n'
        '##JCAMP-DX=5.00\n'
        '##DATA TYPE=CYCLIC VOLTAMMETRY\n'
        + header +
        '##XUNITS=Voltage in V\n'
        '##YUNITS=Current in A\n'
        '##FIRSTX={:.6E}\n'.format(XS[0]) +
        '##LASTX={:.6E}\n'.format(XS[-1]) +
        '##MINX={:.6E}\n'.format(min(XS)) +
        '##MAXX={:.6E}\n'.format(max(XS)) +
        '##MINY={:.6E}\n'.format(min(YS)) +
        '##MAXY={:.6E}\n'.format(max(YS)) +
        '##NPOINTS={}\n'.format(len(XS)) +
        '##FIRSTY={:.6E}\n'.format(YS[0]) +
        '##{}=(XY..XY)\n'.format(table_label) +
        rows +
        '##END=\n'
    )


def _write(tmp_path, name, content):
    path = tmp_path / name
    path.write_text(content)
    return str(path)


def _read_xs(path):
    return np.asarray(JcampNIConverter(JcampBaseConverter(path)).xs)


def test_conformant_xypoints_keeps_real_x(tmp_path):
    """##DATA CLASS=XYPOINTS + ##XYPOINTS= — converter <=1.8.0 and >=1.9.3."""
    path = _write(tmp_path, 'new.jdx', _jcamp('XYPOINTS', 'XYPOINTS'))
    base = JcampBaseConverter(path)

    assert base.data_format == '(XY..XY)'
    np.testing.assert_allclose(_read_xs(path), XS, rtol=1e-6)


def test_mismatched_xydata_label_keeps_real_x(tmp_path):
    """##DATA CLASS=XYPOINTS + ##XYDATA=(XY..XY) — converter 1.9.0-1.9.2.

    Regression guard: before the fix, __set_dataformat could not find a
    'XYPOINTS' key, fell back to '(X++(Y..Y))', and the X axis was
    re-synthesised as linspace(FIRSTX, LASTX) — collapsing a closed cyclic
    sweep to a near-flat line.
    """
    path = _write(tmp_path, 'old.jdx', _jcamp('XYPOINTS', 'XYDATA'))
    base = JcampBaseConverter(path)

    assert base.data_format == '(XY..XY)'
    np.testing.assert_allclose(_read_xs(path), XS, rtol=1e-6)


def test_missing_dataclass_still_finds_xypoints(tmp_path):
    """##XYPOINTS= with no ##DATA CLASS at all."""
    path = _write(tmp_path, 'no_dataclass.jdx', _jcamp(None, 'XYPOINTS'))
    base = JcampBaseConverter(path)

    assert base.dataclass == 'XYPOINTS'
    assert base.data_format == '(XY..XY)'
    np.testing.assert_allclose(_read_xs(path), XS, rtol=1e-6)


def test_equally_spaced_xydata_is_untouched(tmp_path):
    """A genuine ##DATA CLASS=XYDATA + (X++(Y..Y)) file must not be rerouted."""
    content = (
        '##TITLE=Spectrum\n'
        '##JCAMP-DX=5.00\n'
        '##DATA TYPE=CYCLIC VOLTAMMETRY\n'
        '##DATA CLASS=XYDATA\n'
        '##XUNITS=Voltage in V\n'
        '##YUNITS=Current in A\n'
        '##XFACTOR=1.0\n'
        '##YFACTOR=1.0\n'
        '##FIRSTX=0.0\n'
        '##LASTX=6.0\n'
        '##DELTAX=1.0\n'
        '##NPOINTS=7\n'
        '##XYDATA=(X++(Y..Y))\n'
        '0 1 2 3 4 5 6\n'
        '##END=\n'
    )
    path = _write(tmp_path, 'xydata.jdx', content)
    base = JcampBaseConverter(path)

    assert base.dataclass == 'XYDATA_OLD'
    assert base.data_format == '(X++(Y..Y))'
