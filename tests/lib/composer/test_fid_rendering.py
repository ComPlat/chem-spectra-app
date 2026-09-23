"""Bruker FID rendering, which no existing test measures.

`FidBaseConverter` is not fed by data_type.json: it hardcodes typ='NMR',
sets `non_nmr` itself, and carries no spectrum-technique descriptor. So the
composer's descriptor lookup has to fall back on the core's own `non_nmr`
for it. Falling back to the generic curve instead would silently strip the
chemical-shift axis labels off every Bruker render -- and nothing else
would notice: the existing FID tests (test_zip2cv_with_processed_file,
test_zip_to_composer_valid_molfile) assert only how many converters and
composers come back, and never render.

These assert the labels, so that fallback cannot regress unseen.
"""

import pytest
from werkzeug.datastructures import FileStorage

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.model.transformer import TransformerModel
import chem_spectra.lib.composer.technique as technique_module

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
    real_xlabel, real_ylabel = technique_module.plt.xlabel, technique_module.plt.ylabel

    def spy_xlabel(text, **kwargs):
        drawn['x'] = text
        return real_xlabel(text, **kwargs)

    def spy_ylabel(text, **kwargs):
        drawn['y'] = text
        return real_ylabel(text, **kwargs)

    technique_module.plt.xlabel, technique_module.plt.ylabel = spy_xlabel, spy_ylabel
    try:
        composer.tf_img().close()
    finally:
        technique_module.plt.xlabel, technique_module.plt.ylabel = real_xlabel, real_ylabel
    return drawn


def test_fid_core_carries_no_technique_but_reports_nmr(fid_composer):
    """The case the fallback exists for: no descriptor, but non_nmr False."""
    assert getattr(fid_composer.core, 'technique', None) is None
    assert fid_composer.core.non_nmr is False


def test_fid_renders_with_nmr_axis_labels(fid_composer):
    """A descriptor-less core that reports NMR must get the NMR labels.

    These are the values the pre-refactor code produced. A generic-curve
    fallback would give 'X (PPM)' / 'Y (ARBITRARY)' instead.
    """
    drawn = _render_labels(fid_composer)
    assert drawn['x'] == 'Chemical shift (ppm)'
    assert drawn['y'] == 'Intensity (arbitrary)'


def test_technique_fallback_honours_a_cores_own_non_nmr(fid_composer):
    """_technique() resolves a descriptor-less NMR core to the NMR descriptor."""
    assert fid_composer._technique().key == 'NMR'
    assert fid_composer._technique().x_axis == 'chemical_shift'
    assert fid_composer._technique().y_axis == 'intensity'
