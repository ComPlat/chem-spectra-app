from chem_spectra.lib.converter.bagit.base import BagItBaseConverter as BagItConveter
import tempfile
import zipfile
import mimetypes

target_dir = './tests/fixtures/source/bagit/'
cv_layout_path = target_dir + 'cv/File053_BagIt.zip'
aif_layout_path = target_dir + 'aif/aif.zip'
emissions_layout_path = target_dir + 'emissions/emissions.zip'
dls_acf_layout_path = target_dir + 'dls_acf/dls_acf.zip'
dls_intensity_layout_path = target_dir + 'dls_intensity/dls_intensity.zip'
mass_chromatogram_layout_path = target_dir + 'mass_chromatogram/mass_chromatogram.zip'

def assertFileType(file, mimeStr):
    assert mimetypes.guess_type(file.name)[0] == mimeStr

def assertJcampContent(jcamp, field):
    assertFileType(jcamp, 'chemical/x-jcamp-dx')
    jcamp_content = str(jcamp.read())
    assert field in jcamp_content

def test_bagit_converter_failed():
    converter = BagItConveter(None)
    assert converter.data is None
    assert converter.images is None
    assert converter.list_csv is None
    assert converter.combined_image is None

def test_bagit_convert_to_jcamp():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(cv_layout_path, 'r') as z:
            z.extractall(td)

        converter = BagItConveter(td)
        assert converter.data is not None
        assert len(converter.data) == 3

def test_bagit_convert_to_jcamp_cv_layout():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(cv_layout_path, 'r') as z:
            z.extractall(td)

        converter = BagItConveter(td)
        jcamp = converter.data[0]
        assertJcampContent(jcamp, '##DATA TYPE=CYCLIC VOLTAMMETRY')

def test_bagit_convert_to_jcamp_aif_layout():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(aif_layout_path, 'r') as z:
            z.extractall(td)

        converter = BagItConveter(td)
        jcamp = converter.data[0]
        assertJcampContent(jcamp, '##DATA TYPE=SORPTION-DESORPTION MEASUREMENT')

def test_bagit_convert_to_jcamp_emissions_layout():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(emissions_layout_path, 'r') as z:
            z.extractall(td)

        converter = BagItConveter(td)
        jcamp = converter.data[0]
        assertJcampContent(jcamp, '##DATA TYPE=Emissions')

def test_bagit_convert_to_jcamp_dls_acf_layout():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(dls_acf_layout_path, 'r') as z:
            z.extractall(td)

        converter = BagItConveter(td)
        jcamp = converter.data[0]
        assertJcampContent(jcamp, '##DATA TYPE=DLS ACF')

def test_bagit_convert_to_jcamp_dls_intensity_layout():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(dls_intensity_layout_path, 'r') as z:
            z.extractall(td)

        converter = BagItConveter(td)
        jcamp = converter.data[0]
        assertJcampContent(jcamp, '##DATA TYPE=DLS intensity')

def test_bagit_convert_to_jcamp_mass_chromatogram_layout():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(mass_chromatogram_layout_path, 'r') as z:
            z.extractall(td)

        converter = BagItConveter(td)
        jcamp = converter.data[0]
        assertJcampContent(jcamp, '##DATA TYPE=MASS CHROMATOGRAM')

def test_bagit_convert_to_images():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(cv_layout_path, 'r') as z:
            z.extractall(td)

        converter = BagItConveter(td)
        assert converter.images is not None
        assert len(converter.images) == 3
        pngImage = converter.images[0]
        assertFileType(pngImage, 'image/png')
        
def test_bagit_convert_to_csv():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(cv_layout_path, 'r') as z:
            z.extractall(td)

        converter = BagItConveter(td)
        assert converter.list_csv is not None
        assert len(converter.list_csv) == 3
        csvFile = converter.list_csv[0]
        assertFileType(csvFile, 'text/csv')

def test_get_base64_data_failed():
    converter = BagItConveter(None)
    data = converter.get_base64_data()
    assert data is None

def test_get_base64_data_succeed():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(cv_layout_path, 'r') as z:
            z.extractall(td)

        converter = BagItConveter(td)
        list_base64 = converter.get_base64_data()
        assert len(list_base64) == 3
        for base64Str in list_base64:
            assert base64Str.endswith("=")

def test_get_combined_image_failed():
    converter = BagItConveter(None)
    combined_image = converter.combined_image
    assert combined_image is None

def test_get_combined_image():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(cv_layout_path, 'r') as z:
            z.extractall(td)

        converter = BagItConveter(td)
        combined_image = converter.combined_image
        assertFileType(combined_image, 'image/png')

def test_bagit_has_one_file_no_combined_image():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(dls_acf_layout_path, 'r') as z:
            z.extractall(td)

        converter = BagItConveter(td)
        assert converter.combined_image is None
