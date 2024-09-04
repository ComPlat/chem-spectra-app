import tempfile
import zipfile
import mimetypes
import base64
from chem_spectra.lib.converter.lcms.base import LCMSBaseConverter as LCMSConveter

target_dir = './tests/fixtures/source/lcms/lcms.zip'

def assertFileType(file, mimeStr):
    assert mimetypes.guess_type(file.name)[0] == mimeStr

def assertJcampContent(jcamp, field):
    assertFileType(jcamp, 'chemical/x-jcamp-dx')
    jcamp_content = str(jcamp.read())
    assert field in jcamp_content
    
def isBase64(encodeString):
    plainStr = base64.b64decode(encodeString)
    encodedStr = base64.b64encode(plainStr).decode("utf-8")
    assert encodedStr == encodeString

def test_lcms_converter_failed():
    converter = LCMSConveter(None)
    assert converter.data is None

def test_lcms_converter_success():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(target_dir, 'r') as z:
            z.extractall(td)

        converter = LCMSConveter(td, fname='lcms')
        assert converter.data is not None
        assert len(converter.data) == 5

def test_lcms_converter_tic_positive():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(target_dir, 'r') as z:
            z.extractall(td)

        converter = LCMSConveter(td, fname='lcms')
        tic_positive = converter.data[0]
        assert len(tic_positive['time']) > 0
        assert len(tic_positive['Intensity']) > 0

def test_lcms_converter_tic_negative():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(target_dir, 'r') as z:
            z.extractall(td)

        converter = LCMSConveter(td, fname='lcms')
        tic_positive = converter.data[1]
        assert len(tic_positive['time']) > 0
        assert len(tic_positive['Intensity']) > 0

def test_lcms_converter_uvvis():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(target_dir, 'r') as z:
            z.extractall(td)

        converter = LCMSConveter(td, fname='lcms')
        uvvis = converter.data[2]
        for wavelength in uvvis.keys():
            assert len(uvvis[wavelength]['RetentionTime']) > 0
            assert len(uvvis[wavelength]['DetectorSignal']) > 0

def test_lcms_converter_spectra_positive():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(target_dir, 'r') as z:
            z.extractall(td)

        converter = LCMSConveter(td, fname='lcms')
        spectra_positive = converter.data[3]
        for time in spectra_positive.keys():
            assert len(spectra_positive[time]['mz']) > 0
            assert len(spectra_positive[time]['intensities']) > 0

def test_lcms_converter_spectra_negative():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(target_dir, 'r') as z:
            z.extractall(td)

        converter = LCMSConveter(td, fname='lcms')
        spectra_negative = converter.data[4]
        for time in spectra_negative.keys():
            assert len(spectra_negative[time]['mz']) > 0
            assert len(spectra_negative[time]['intensities']) > 0

