import pytest
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter

source_nmr = './tests/fixtures/source/1H.dx'

@pytest.fixture
def jcamp_file_1h():
    return source_nmr

def test_init_jcamp_ni_converter_failed():
    with pytest.raises(Exception) as error:
        _ = JcampNIConverter(None)
        
    assert error is not None

def test_init_jcamp_ni_success(jcamp_file_1h):
    base_converter = JcampBaseConverter(jcamp_file_1h)
    ni_converter = JcampNIConverter(base=base_converter)
    
    assert ni_converter is not None
    assert ni_converter.base == base_converter
    