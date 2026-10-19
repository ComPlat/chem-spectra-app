"""A FID must reach JcampTechniqueConverter carrying its descriptor.

`FidBaseConverter` hardcodes typ='NMR' and is not fed by data_type.json, so
nothing derives a descriptor for it. But it is handed straight to
JcampTechniqueConverter (model/transformer.py:213), which copies
`getattr(base, 'technique', None)` across -- so without one set here, every
Bruker render travelled with `technique = None`.

That was survivable only because the composer falls back on the core's own
`non_nmr`. Any consumer reading a field off the descriptor -- the threshold
is the first -- would raise AttributeError on the FID path instead.
"""

import pytest
from werkzeug.datastructures import FileStorage

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.lib.converter.jcamp.techniques import SPECTRUM_TECHNIQUES
from chem_spectra.model.transformer import TransformerModel

source_zip = './tests/fixtures/source/bruker/1H.zip'
source_molfile = './tests/fixtures/source/molfile/svs813f1_B.mol'


@pytest.fixture
def fid_converter():
    with open(source_zip, 'rb') as handle:
        container = FileContainer(FileStorage(handle))
        with open(source_molfile, 'rb') as molhandle:
            molfile = FileContainer(FileStorage(molhandle))
            converters, _, _ = TransformerModel(
                container, molfile=molfile, params={'ext': 'zip'},
            ).zip2cvp()
    return converters[0] if isinstance(converters, list) else converters


def test_fid_converter_carries_the_nmr_descriptor(fid_converter):
    assert fid_converter.technique is SPECTRUM_TECHNIQUES['NMR']


def test_fid_descriptor_survives_into_the_technique_converter(fid_converter):
    """The copy-across in JcampTechniqueConverter.__init__ must not yield None."""
    assert fid_converter.technique is not None
    assert fid_converter.technique.key == 'NMR'
