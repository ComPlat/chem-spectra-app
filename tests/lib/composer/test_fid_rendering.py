"""Bruker FID rendering, which no existing test measures.

`FidBaseConverter` is not fed by data_type.json: it hardcodes typ='NMR' and
sets the NMR descriptor itself. Resolving to the generic curve instead
would silently strip the chemical-shift axis labels off every Bruker render
-- and nothing else would notice: the existing FID tests
(test_zip2cv_with_processed_file, test_zip_to_composer_valid_molfile) assert
only how many converters and composers come back, and never render.

These assert the labels, so that cannot regress unseen.
"""

import pytest
from werkzeug.datastructures import FileStorage

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.model.transformer import TransformerModel
import chem_spectra.lib.composer.curve as curve_module

source_zip = './tests/fixtures/source/bruker/1H.zip'
source_molfile = './tests/fixtures/source/molfile/svs813f1_B.mol'


@pytest.fixture
def fid_composer():
    with open(source_zip, 'rb') as handle:
        container = FileContainer(FileStorage(handle))
        with open(source_molfile, 'rb') as molhandle:
            molfile = FileContainer(FileStorage(molhandle))
            _, composers, _ = TransformerModel(
                container, molfile=molfile, params={'ext': 'zip'},
            ).zip2cvp()
    return composers[0] if isinstance(composers, list) else composers


def _render_labels(composer):
    drawn = {}
    real_xlabel, real_ylabel = curve_module.plt.xlabel, curve_module.plt.ylabel

    def spy_xlabel(text, **kwargs):
        drawn['x'] = text
        return real_xlabel(text, **kwargs)

    def spy_ylabel(text, **kwargs):
        drawn['y'] = text
        return real_ylabel(text, **kwargs)

    curve_module.plt.xlabel, curve_module.plt.ylabel = spy_xlabel, spy_ylabel
    try:
        composer.tf_img().close()
    finally:
        curve_module.plt.xlabel, curve_module.plt.ylabel = real_xlabel, real_ylabel
    return drawn


def test_fid_core_carries_the_nmr_descriptor(fid_composer):
    """FidBaseConverter sets the descriptor itself; typ is 'NMR' by construction.

    This used to assert the opposite -- that the core carried no descriptor
    at all -- which pinned an accidental absence rather than behaviour. Every
    core now carries one (tests/lib/converter/
    test_every_core_carries_a_descriptor.py).
    """
    assert fid_composer.core.technique is not None
    assert fid_composer.core.technique.key == 'NMR'


def test_fid_renders_with_nmr_axis_labels(fid_composer):
    """A FID gets the NMR labels.

    These are the values the pre-refactor code produced. A generic-curve
    fallback would give 'X (PPM)' / 'Y (ARBITRARY)' instead.
    """
    drawn = _render_labels(fid_composer)
    assert drawn['x'] == 'Chemical shift (ppm)'
    assert drawn['y'] == 'Intensity (arbitrary)'


def test_a_fid_resolves_to_the_nmr_descriptor(fid_composer):
    assert fid_composer._technique().key == 'NMR'
    assert fid_composer._technique().x_axis == 'chemical_shift'
    assert fid_composer._technique().y_axis == 'intensity'


def test_a_core_without_a_descriptor_is_not_taken_for_nmr(fid_composer):
    """`non_nmr = False` no longer stands in for a descriptor.

    It used to: `_technique()` resolved such a core to NMR, and the NMRium
    converter relied on that. Every core now carries its descriptor
    (tests/lib/converter/test_every_core_carries_a_descriptor.py), so a core
    without one is the generic curve, whatever flags it sets.
    """
    class DescriptorlessCore:
        non_nmr = False

    composer = object.__new__(type(fid_composer))
    composer.core = DescriptorlessCore()
    assert composer._technique().key == ''
