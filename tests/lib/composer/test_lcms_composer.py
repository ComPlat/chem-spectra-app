import json
import pytest
import tempfile
import zipfile
from chem_spectra.lib.converter.lcms.base import LCMSBaseConverter
from chem_spectra.lib.composer.lcms import LCMSComposer

source = './tests/fixtures/source/lcms/lcms.zip'

@pytest.fixture
def zip_file():
    return source

def test_init_lcms_composer_failed():
    with pytest.raises(Exception) as error:
        _ = LCMSComposer(None)
        
    assert error is not None

def test_init_lcms_composer_success(zip_file):
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(zip_file, 'r') as z:
            z.extractall(td)
          
        lcms_converter = LCMSBaseConverter(td)
        lcms_composer = LCMSComposer(core=lcms_converter)
        
        assert lcms_composer is not None
        assert lcms_composer.core == lcms_converter

# def test_ms_composer_original_metadata(jcamp_file):
#     base_converter = JcampBaseConverter(jcamp_file)
#     ms_converter = JcampMSConverter(base=base_converter)
#     ms_composer = MSComposer(core=ms_converter)

#     assert ms_composer is not None
#     assert '$$ === CHEMSPECTRA ORIGINAL METADATA ===\n' in ms_composer.meta
