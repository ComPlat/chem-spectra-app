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
    

def test_missing_datatype_header_does_not_raise(tmp_path):
    """A JCAMP with no ##DATA TYPE= header at all.

    `self.dic['DATATYPE']` was unguarded, so such a file raised
    KeyError('DATATYPE') straight out of the request.
    """
    body = '\n'.join(
        line for line in open(source_ir, encoding='utf-8').read().split('\n')
        if 'DATA TYPE' not in line
    )
    target = tmp_path / 'no_datatype.jdx'
    target.write_text(body, encoding='utf-8')

    converter = JcampBaseConverter(str(target))
    assert converter.datatypes == []
    assert converter.non_nmr is False   # master's behaviour for typ == ''
