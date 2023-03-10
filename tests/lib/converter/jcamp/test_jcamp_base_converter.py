import pytest
from werkzeug.datastructures import FileStorage
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter

source_nmr = './tests/fixtures/source/1H.dx'
source_ir = './tests/fixtures/source/IR.dx'

@pytest.fixture
def jcamp_file_1h():
    # with open(source, 'rb') as f:
    #     jcamp = FileContainer(FileStorage(f))
    # return jcamp
    return source_nmr

@pytest.fixture
def jcamp_file_ir():
    return source_ir

def test_init_jcamp_base_converter_failed():
    with pytest.raises(Exception) as error:
        _ = JcampBaseConverter(path="ddd")
        
    assert error is not None

def test_init_jcamp_base_converter_success(jcamp_file_1h):
    converter = JcampBaseConverter(jcamp_file_1h)
    assert converter is not None

def test_convert_jcamp_nmr(jcamp_file_1h):
    converter = JcampBaseConverter(jcamp_file_1h)
    assert converter is not None
    assert converter.non_nmr == False
    assert converter.ncl == "1H"
  
def test_convert_jcamp_non_nmr(jcamp_file_ir):
    converter = JcampBaseConverter(jcamp_file_ir)
    assert converter is not None
    assert converter.non_nmr == True
    