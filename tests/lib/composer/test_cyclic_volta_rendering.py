"""Cyclic voltammetry rendering, pinned before the dispatch refactor.

`is_cyclic_volta` is the heaviest of the sixteen flags (22 references) and
carries the most behaviour: y-scaling from `yScaleFactor`, a scientific
axis exponent derived from the data, raw axis labels rather than the
generic "X (unit)" form, and max/min peak pairs with `isRef` handling.

The existing tests/test_cyclic_volta.py covers the flag, the parsed x/y and
one metadata string. None of the rendering behaviour was covered, which is
exactly what a registry field named `cv_scaling` would have to reproduce.
"""

import json

import numpy as np
import pytest

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter
from chem_spectra.lib.composer.ni import NIComposer
import chem_spectra.lib.composer.ni as ni_module

SOURCE = './tests/fixtures/source/cyclicvoltammetry/RCV_LSH-R444_full+Fc.jdx'


def _composer(params=False):
    return NIComposer(JcampNIConverter(JcampBaseConverter(SOURCE, params)))


def _render_capture(composer):
    """Run tf_img, capturing the first plotted series and the axis labels."""
    drawn = {}
    real_plot = ni_module.plt.plot
    real_xlabel, real_ylabel = ni_module.plt.xlabel, ni_module.plt.ylabel

    def spy_plot(*args, **kwargs):
        if 'series' not in drawn and len(args) >= 2:
            drawn['series'] = np.asarray(args[1], dtype=float)
        return real_plot(*args, **kwargs)

    def spy_xlabel(text, **kwargs):
        drawn['xlabel'] = text
        return real_xlabel(text, **kwargs)

    def spy_ylabel(text, **kwargs):
        drawn['ylabel'] = text
        return real_ylabel(text, **kwargs)

    ni_module.plt.plot = spy_plot
    ni_module.plt.xlabel = spy_xlabel
    ni_module.plt.ylabel = spy_ylabel
    try:
        composer.tf_img().close()
    finally:
        ni_module.plt.plot = real_plot
        ni_module.plt.xlabel = real_xlabel
        ni_module.plt.ylabel = real_ylabel
    return drawn


def _cv_params(scale):
    # the caller-facing field is `cyclic_volta`; parse_params json-loads it and
    # republishes it as `cyclicvolta`, and it dereferences spectraList
    # unconditionally, so that key has to be present
    return {'cyclic_volta': json.dumps({
        'spectraList': [{'list': [], 'shift': {}, 'hasRefPeak': False}],
        'cvDisplay': {'yScaleFactor': scale},
    })}


def test_cv_labels_are_raw_not_wrapped():
    """CV is the only technique whose axis labels are printed verbatim.

    Every other non-NMR technique gets "X (unit)" / "Y (unit)".
    """
    drawn = _render_capture(_composer())
    assert not drawn['xlabel'].startswith('X (')
    assert not drawn['ylabel'].startswith('Y (')
    assert drawn['xlabel'] == 'V vs Ref'   # the raw ##XUNITS value


def test_cv_without_scaling_plots_the_raw_series():
    composer = _composer()
    drawn = _render_capture(composer)
    assert composer._cv_density_scale == 1.0
    np.testing.assert_allclose(drawn['series'][:5], composer.core.ys[:5])


@pytest.mark.parametrize('scale', [2.0, 0.5])
def test_cv_yscalefactor_scales_the_plotted_series(scale):
    composer = _composer(_cv_params(scale))
    raw = np.asarray(composer.core.ys, dtype=float).copy()
    drawn = _render_capture(composer)

    assert composer._cv_density_scale == scale
    np.testing.assert_allclose(drawn['series'][:5], raw[:5] * scale)
    # the underlying data is untouched; only the drawn series is scaled
    np.testing.assert_allclose(np.asarray(composer.core.ys[:5], dtype=float), raw[:5])


def test_cv_axis_exponent_follows_the_data_magnitude():
    """The y axis is drawn as a mantissa plus a shared 10^n label."""
    composer = _composer()
    _render_capture(composer)
    expected = int(np.floor(np.log10(np.max(np.abs(composer.core.ys)))))
    assert composer._cv_axis_exp == expected
    assert composer._cv_axis_base == 10.0 ** expected


def test_cv_malformed_yscalefactor_falls_back_to_one():
    composer = _composer({'cyclic_volta': json.dumps({
        'spectraList': [{'list': [], 'shift': {}, 'hasRefPeak': False}],
        'cvDisplay': {'yScaleFactor': 'not-a-number'},
    })})
    _render_capture(composer)
    assert composer._cv_density_scale == 1.0
