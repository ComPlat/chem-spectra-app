"""Bruker FID rendering, which no existing test measures.

`FidBaseConverter` is not fed by data_type.json: it hardcodes typ='NMR' and
sets `non_nmr` itself. Everything about how its output is drawn therefore
rides on that one flag.

The existing FID tests (test_zip2cv_with_processed_file,
test_zip_to_composer_valid_molfile) assert only how many converters and
composers come back and never render, so nothing measures the axis labels.
These do, before any dispatch logic is moved.
"""

import pytest
from werkzeug.datastructures import FileStorage

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.model.transformer import TransformerModel
import chem_spectra.lib.composer.ni as ni_module

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
    real_xlabel, real_ylabel = ni_module.plt.xlabel, ni_module.plt.ylabel

    def spy_xlabel(text, **kwargs):
        drawn['x'] = text
        return real_xlabel(text, **kwargs)

    def spy_ylabel(text, **kwargs):
        drawn['y'] = text
        return real_ylabel(text, **kwargs)

    ni_module.plt.xlabel, ni_module.plt.ylabel = spy_xlabel, spy_ylabel
    try:
        composer.tf_img().close()
    finally:
        ni_module.plt.xlabel, ni_module.plt.ylabel = real_xlabel, real_ylabel
    return drawn



def test_fid_renders_with_nmr_axis_labels(fid_composer):
    """A core that reports non_nmr False must get the NMR axis labels."""
    drawn = _render_labels(fid_composer)
    assert drawn['x'] == 'Chemical shift (ppm)'
    assert drawn['y'] == 'Intensity (arbitrary)'




def test_fid_core_reports_nmr():
    """The single flag everything downstream keys off."""
    with open(source_zip, 'rb') as handle:
        container = FileContainer(FileStorage(handle))
        with open(source_molfile, 'rb') as molhandle:
            molfile = FileContainer(FileStorage(molhandle))
            _, composers, _ = TransformerModel(
                container, molfile=molfile, params={'ext': 'zip'},
            ).zip2cvp()
    composer = composers[0] if isinstance(composers, list) else composers
    assert composer.core.non_nmr is False
