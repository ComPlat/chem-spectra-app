"""Coverage for TransformerModel.to_converter().

It had none. `cv, _ = self.jcamp2cvp()` unpacked two values from a
three-value return, so every JCAMP file raised ValueError on the only
caller, /predict/by_peaks_form (controller/inference_api.py).
"""

import pytest
from werkzeug.datastructures import FileStorage

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter
from chem_spectra.lib.converter.jcamp.ms import JcampMSConverter
from chem_spectra.model.transformer import TransformerModel

source_nmr = './tests/fixtures/source/1H.dx'
source_ir = './tests/fixtures/source/IR.dx'
source_ms = './tests/fixtures/source/ms/MS_ESI.jdx'


def _converter_for(path, params):
    with open(path, 'rb') as handle:
        file_container = FileContainer(FileStorage(handle))
        return TransformerModel(
            file_container, molfile=None, params=params,
        ).to_converter()


@pytest.mark.parametrize('path', [source_nmr, source_ir])
def test_to_converter_returns_a_converter_for_jcamp(path):
    # this raised ValueError: too many values to unpack (expected 2)
    assert isinstance(_converter_for(path, {'ext': 'jdx'}), JcampTechniqueConverter)


def test_to_converter_returns_ms_converter_for_ms_jcamp():
    assert isinstance(_converter_for(source_ms, {'ext': 'jdx'}), JcampMSConverter)


def test_to_converter_and_to_composer_agree_on_the_jcamp_path():
    """Both dispatch on the same extension rules and must not diverge.

    The duplicated dispatch in these two methods is what hid the unpack bug:
    to_composer unpacked three values correctly, to_converter did not.
    """
    with open(source_nmr, 'rb') as handle:
        file_container = FileContainer(FileStorage(handle))
        model = TransformerModel(file_container, molfile=None, params={'ext': 'jdx'})
        converter = model.to_converter()

    with open(source_nmr, 'rb') as handle:
        file_container = FileContainer(FileStorage(handle))
        model = TransformerModel(file_container, molfile=None, params={'ext': 'jdx'})
        composer, _ = model.to_composer()

    assert converter.typ == composer.core.typ
