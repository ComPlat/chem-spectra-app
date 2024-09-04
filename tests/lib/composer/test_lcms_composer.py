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
        assert len(lcms_composer.data) == 5

def test_lcms_composer_tic_postive(zip_file):
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(zip_file, 'r') as z:
            z.extractall(td)
          
        lcms_converter = LCMSBaseConverter(td)
        lcms_composer = LCMSComposer(core=lcms_converter)
        
        tic_postive_jcamp = lcms_composer.data[0]
        assert tic_postive_jcamp is not None
        with open(tic_postive_jcamp.name) as file:
            file_content = file.read()
            assert file_content != ""
            assert '##$CSCATEGORY=TIC SPECTRUM' in file_content

def test_lcms_composer_tic_negative(zip_file):
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(zip_file, 'r') as z:
            z.extractall(td)
          
        lcms_converter = LCMSBaseConverter(td)
        lcms_composer = LCMSComposer(core=lcms_converter)
        
        tic_postive_jcamp = lcms_composer.data[1]
        assert tic_postive_jcamp is not None
        with open(tic_postive_jcamp.name) as file:
            file_content = file.read()
            assert file_content != ""
            assert '##$CSCATEGORY=TIC SPECTRUM' in file_content

def test_lcms_composer_tic_uvvis(zip_file):
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(zip_file, 'r') as z:
            z.extractall(td)
          
        lcms_converter = LCMSBaseConverter(td)
        lcms_composer = LCMSComposer(core=lcms_converter)
        
        uvvis_jcamp = lcms_composer.data[2]
        assert uvvis_jcamp is not None
        with open(uvvis_jcamp.name) as file:
            file_content = file.read()
            assert file_content != ""
            assert '##$CSCATEGORY=UVVIS SPECTRUM' in file_content
