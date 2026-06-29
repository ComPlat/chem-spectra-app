import io
import tempfile
import zipfile
import os

import pytest
from werkzeug.datastructures import FileStorage
from chem_spectra.lib.converter.fid.base import FidBaseConverter
from chem_spectra.lib.converter.fid.bruker import FidHasBruckerProcessed
from chem_spectra.lib.converter.bagit.base import BagItBaseConverter
from chem_spectra.model.transformer import TransformerModel
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter
from chem_spectra.lib.composer.ni import NIComposer

source_dir_1h_bruker = './tests/fixtures/source/bruker/1H.zip'
source_dir_1h_jcamp = './tests/fixtures/source/1H.dx'
filename_1h = '1H'
params_1h_bruker = {'fname':'1H.zip', 'ext':'zip'}
params_1h_jcamp = {'fname':'1H.dx', 'ext':'dx'}
params_1h_bruker_check_nmr = {'fname':'1H.zip', 'ext':'zip', 'simulatenmr': True}
params_1h_jcamp_check_nmr = {'fname':'1H.dx', 'ext':'dx', 'simulatenmr': True}

source_dir_molfile = './tests/fixtures/source/molfile/svs813f1_B.mol'
source_dir_invalid_molfile = './tests/fixtures/source/molfile/invalid_molfile.mol'

def is_list_of_instance(list_data, cls):
    for item in list_data:
        if (isinstance(item, cls) == False):
            return False
    return True

def test_init_with_single_file():
    file = open(source_dir_1h_jcamp, "rb")
    molfile = open(source_dir_molfile, "r")
    tranform_model = TransformerModel(file, molfile=molfile, params=params_1h_jcamp)
    
    assert tranform_model != None
    
def test_init_with_multiple_files():
    file_1 = open(source_dir_1h_jcamp, "rb")
    file_2 = open(source_dir_1h_jcamp, "rb")
    molfile = open(source_dir_molfile, "r")
    tranform_model = TransformerModel(None, molfile=molfile, params=params_1h_jcamp, multiple_files=[file_1, file_2])
    
    assert tranform_model != None
    
def test_combine():
    with open(source_dir_1h_jcamp, 'rb') as f:
        file_1 = FileContainer(FileStorage(f))

    with open(source_dir_1h_jcamp, 'rb') as f:
        file_2 = FileContainer(FileStorage(f))

    molfile = open(source_dir_molfile, "r")
    tranform_model = TransformerModel(None, molfile=molfile, params=params_1h_jcamp, multiple_files=[file_1, file_2])
    
    assert tranform_model is not None
    
    tf = tranform_model.tf_combine()
    assert tf is not None


def test_zip2cv_with_processed_file():
    with open(source_dir_molfile, 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    with open(source_dir_1h_bruker, 'rb') as f:
        file = FileContainer(FileStorage(f))

    tranform_model = TransformerModel(file, molfile=molfile, params=params_1h_bruker)

    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(source_dir_1h_bruker, 'r') as z:
            z.extractall(td)
        
        target_dir = os.path.join(td, '1')

        list_converters, list_composers, invalid_molfile = tranform_model.zip2cv_with_processed_file(target_dir=target_dir, params=params_1h_bruker, file_name=filename_1h)
        assert len(list_converters) == 2
        assert len(list_composers) == 2
        assert invalid_molfile is False

def test_zip_to_composer_invalid_molfile():
    with open(source_dir_invalid_molfile, 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    with open(source_dir_1h_bruker, 'rb') as f:
        file = FileContainer(FileStorage(f))
    
    tranform_model = TransformerModel(file, molfile=molfile, params=params_1h_bruker_check_nmr)

    list_converters, list_composers, invalid_molfile = tranform_model.zip2cvp()
    assert len(list_converters) == 2
    assert len(list_composers) == 2
    assert invalid_molfile is True

def test_zip_to_composer_valid_molfile():
    with open(source_dir_molfile, 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    with open(source_dir_1h_bruker, 'rb') as f:
        file = FileContainer(FileStorage(f))
    
    tranform_model = TransformerModel(file, molfile=molfile, params=params_1h_bruker)

    list_converters, list_composers, invalid_molfile = tranform_model.zip2cvp()
    assert len(list_converters) == 2
    assert len(list_composers) == 2
    assert invalid_molfile is False

def test_jcamp_to_composer_invalid_molfile():
    with open(source_dir_invalid_molfile, 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    with open(source_dir_1h_jcamp, 'rb') as f:
        file = FileContainer(FileStorage(f))

    tranform_model = TransformerModel(file, molfile=molfile, params=params_1h_jcamp_check_nmr)
    converter, composer, invalid_molfile = tranform_model.jcamp2cvp()
    assert isinstance(converter, JcampNIConverter)
    assert isinstance(composer, NIComposer)
    assert invalid_molfile is True

def test_jcamp_to_composer_valid_molfile():
    with open(source_dir_molfile, 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    with open(source_dir_1h_jcamp, 'rb') as f:
        file = FileContainer(FileStorage(f))

    tranform_model = TransformerModel(file, molfile=molfile, params=params_1h_jcamp)
    converter, composer, invalid_molfile = tranform_model.jcamp2cvp()
    assert isinstance(converter, JcampNIConverter)
    assert isinstance(composer, NIComposer)
    assert invalid_molfile is False

def test_to_composer_jcamp_invalid_molfile():
    with open(source_dir_invalid_molfile, 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    with open(source_dir_1h_jcamp, 'rb') as f:
        file = FileContainer(FileStorage(f))
    
    tranform_model = TransformerModel(file, molfile=molfile, params=params_1h_jcamp_check_nmr)

    composer, invalid_molfile = tranform_model.to_composer()
    assert isinstance(composer, NIComposer)
    assert invalid_molfile is True


def test_to_composer_zip_invalid_molfile():
    with open(source_dir_invalid_molfile, 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    with open(source_dir_1h_bruker, 'rb') as f:
        file = FileContainer(FileStorage(f))
    
    tranform_model = TransformerModel(file, molfile=molfile, params=params_1h_bruker_check_nmr)

    list_composers, invalid_molfile = tranform_model.to_composer()
    assert is_list_of_instance(list_composers, NIComposer)
    assert invalid_molfile is True

def test_to_composer_jcamp_valid_molfile():
    with open(source_dir_molfile, 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    with open(source_dir_1h_jcamp, 'rb') as f:
        file = FileContainer(FileStorage(f))
    
    tranform_model = TransformerModel(file, molfile=molfile, params=params_1h_jcamp)

    composer, invalid_molfile = tranform_model.to_composer()
    assert isinstance(composer, NIComposer)
    assert invalid_molfile is False

def test_to_composer_zip_valid_molfile():
    with open(source_dir_molfile, 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    with open(source_dir_1h_bruker, 'rb') as f:
        file = FileContainer(FileStorage(f))
    
    tranform_model = TransformerModel(file, molfile=molfile, params=params_1h_bruker)

    list_composers, invalid_molfile = tranform_model.to_composer()
    assert is_list_of_instance(list_composers, NIComposer)
    assert invalid_molfile is False

def test_to_composer_mass_spectrum():
    # TODO: implement later
    pass

def test_to_composer_cdf():
    # TODO: implement later
    pass

def test_convert2jcamp():
    with open(source_dir_molfile, 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    with open(source_dir_1h_jcamp, 'rb') as f:
        file = FileContainer(FileStorage(f))
    
    tranform_model = TransformerModel(file, molfile=molfile, params=params_1h_jcamp)

    jcamp = tranform_model.convert2jcamp()
    assert jcamp is not None

def test_convert2jcamp_bagit():
    # TODO: implement later
    pass

def test_convert2img():
    with open(source_dir_molfile, 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    with open(source_dir_1h_jcamp, 'rb') as f:
        file = FileContainer(FileStorage(f))
    
    tranform_model = TransformerModel(file, molfile=molfile, params=params_1h_jcamp)

    image = tranform_model.convert2img()
    assert image is not None

def test_convert2img_bagit():
    # TODO: implement later
    pass

def test_convert2jcamp_img():
    with open(source_dir_molfile, 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    with open(source_dir_1h_jcamp, 'rb') as f:
        file = FileContainer(FileStorage(f))
    
    tranform_model = TransformerModel(file, molfile=molfile, params=params_1h_jcamp)

    jcamp, image, csv = tranform_model.convert2jcamp_img()
    assert jcamp is not None
    assert image is not None
    assert csv is None

def test_convert2jcamp_img_bagit():
    # TODO: implement later
    pass


_LCMS_FLAT_JDX = """##TITLE={title}
##JCAMP-DX=5.00 $$ chemotion-converter-app (0.1.0)
##DATA TYPE=LC/MS
##DATA CLASS=XYPOINTS
##XUNITS=Time (s)
##YUNITS=Intensity
##NPOINTS=3
##XYPOINTS=(XY..XY)
1.0, 100.0
2.0, 200.0
3.0, 300.0
##END=
"""


def _make_zip_filecontainer(name, files):
    """Return a FileContainer wrapping an in-memory zip with ``files``.

    ``files`` is a mapping ``arcname → text content``.
    """
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED) as z:
        for arcname, content in files.items():
            z.writestr(arcname, content)
    buffer.seek(0)
    storage = FileStorage(stream=buffer, filename=name, content_type='application/zip')
    return FileContainer(storage)


def test_zip2cvp_with_flat_lcms_jdx_files_routes_to_bagit():
    """A plain zip containing LCMS ``.jdx`` files at the root must be picked
    up and routed through ``BagItBaseConverter`` (no BagIt envelope needed)."""
    file = _make_zip_filecontainer('lcms.zip', {
        'tic.peak.jdx': _LCMS_FLAT_JDX.format(title='TIC peak'),
        'tic.edit.jdx': _LCMS_FLAT_JDX.format(title='TIC edit'),
        'ms.peak.jdx': _LCMS_FLAT_JDX.format(title='MS peak'),
        'ms.edit.jdx': _LCMS_FLAT_JDX.format(title='MS edit'),
    })
    transformer = TransformerModel(file, params={'fname': 'lcms.zip', 'ext': 'zip'})

    converter, composer, invalid_molfile = transformer.zip2cvp()

    assert isinstance(converter, BagItBaseConverter)
    assert composer is converter
    assert invalid_molfile is False
    assert converter.data is not None
    assert len(converter.data) == 4


def test_zip2cvp_with_flat_lcms_in_single_subfolder():
    """LCMS exports often unzip into a single named subfolder. We must look
    one level deep when the root has no .jdx."""
    file = _make_zip_filecontainer('lcms_subdir.zip', {
        'export_2024/tic.peak.jdx': _LCMS_FLAT_JDX.format(title='TIC peak'),
        'export_2024/ms.peak.jdx': _LCMS_FLAT_JDX.format(title='MS peak'),
    })
    transformer = TransformerModel(file, params={'fname': 'lcms_subdir.zip', 'ext': 'zip'})

    converter, _, _ = transformer.zip2cvp()

    assert isinstance(converter, BagItBaseConverter)
    assert converter.data is not None
    assert len(converter.data) == 2


def test_zip2cvp_returns_false_when_no_recognised_content():
    """A zip with neither Bruker FID, BagIt envelope, nor any .jdx must
    keep returning the historical (False, False, False) tuple so existing
    callers behave unchanged."""
    file = _make_zip_filecontainer('garbage.zip', {
        'README.txt': 'not a spectrum',
        'notes.md': '# nothing here',
    })
    transformer = TransformerModel(file, params={'fname': 'garbage.zip', 'ext': 'zip'})

    converter, composer, invalid_molfile = transformer.zip2cvp()

    assert converter is False
    assert composer is False
    assert invalid_molfile is False


def test_zip2cvp_returns_false_on_zip_slip_attempt():
    """A zip file containing directory traversal paths (Zip Slip) must fail
    closed by returning (False, False, False) rather than raising an uncaught exception."""
    file = _make_zip_filecontainer('zip_slip.zip', {
        '../../escape.txt': 'malicious content',
    })
    transformer = TransformerModel(file, params={'fname': 'zip_slip.zip', 'ext': 'zip'})
    converter, composer, invalid_molfile = transformer.zip2cvp()
    assert converter is False
    assert composer is False
    assert invalid_molfile is False
