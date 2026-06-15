import pytest
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter

from tests.dataset_catalog import dataset_path_str

source_nmr = dataset_path_str('NMR-021')

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

def test_init_jcamp_ni_nmr_label(jcamp_file_1h):
    base_converter = JcampBaseConverter(jcamp_file_1h)
    ni_converter = JcampNIConverter(base=base_converter)

    assert ni_converter.label == {'x': 'PPM', 'y': 'ARBITRARY'}
    