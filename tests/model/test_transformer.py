import tempfile
import zipfile
import os

import pytest
from werkzeug.datastructures import FileStorage
from chem_spectra.lib.converter.fid.base import FidBaseConverter
from chem_spectra.lib.converter.fid.bruker import FidHasBruckerProcessed
from chem_spectra.model.transformer import TransformerModel
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter
from chem_spectra.lib.composer.ni import NIComposer

source_dir_1h_bruker = './tests/fixtures/source/bruker/1H.zip'
source_dir_1h_jcamp = './tests/fixtures/source/1H.dx'
filename_1h = '1H'
params_1h_bruker = {'fname':'1H.zip', 'ext':'zip'}
params_1h_jcamp = {'fname':'1H.dx', 'ext':'dx'}

source_dir_molfile = './tests/fixtures/source/molfile/svs813f1_B.mol'
source_dir_invalid_molfile = './tests/fixtures/source/molfile/invalid_molfile.mol'

def is_list_of_instance(list_data, cls):
    for item in list_data:
        if (isinstance(item, cls) == False):
            return False
    return True

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
    
    tranform_model = TransformerModel(file, molfile=molfile, params=params_1h_bruker)

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

    tranform_model = TransformerModel(file, molfile=molfile, params=params_1h_jcamp)
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
    
    tranform_model = TransformerModel(file, molfile=molfile, params=params_1h_jcamp)

    composer, invalid_molfile = tranform_model.to_composer()
    assert isinstance(composer, NIComposer)
    assert invalid_molfile is True


def test_to_composer_zip_invalid_molfile():
    with open(source_dir_invalid_molfile, 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    with open(source_dir_1h_bruker, 'rb') as f:
        file = FileContainer(FileStorage(f))
    
    tranform_model = TransformerModel(file, molfile=molfile, params=params_1h_bruker)

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
