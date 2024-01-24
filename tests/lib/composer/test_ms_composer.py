import json
import pytest
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ms import JcampMSConverter
from chem_spectra.lib.composer.ms import MSComposer

source = './tests/fixtures/source/ms/ms_v6.dx'

@pytest.fixture
def jcamp_file():
    return source

def test_init_ms_composer_failed():
    with pytest.raises(Exception) as error:
        _ = MSComposer(None)
        
    assert error is not None

def test_init_ms_composer_success(jcamp_file):
    base_converter = JcampBaseConverter(jcamp_file)
    ms_converter = JcampMSConverter(base=base_converter)
    ms_composer = MSComposer(core=ms_converter)
    
    assert ms_composer is not None
    assert ms_composer.core == ms_converter

def test_ms_composer_original_metadata(jcamp_file):
    base_converter = JcampBaseConverter(jcamp_file)
    ms_converter = JcampMSConverter(base=base_converter)
    ms_composer = MSComposer(core=ms_converter)

    assert ms_composer is not None
    assert '$$ === CHEMSPECTRA ORIGINAL METADATA ===\n' in ms_composer.meta
