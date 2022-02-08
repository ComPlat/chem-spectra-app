from werkzeug.datastructures import FileStorage
import tempfile
import zipfile
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.lib.shared.buffer import store_str_in_tmp, store_byte_in_tmp

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.bagit.base import BagItBaseConverter

target_dir = './tests/fixtures/'
source_dir = 'source/cyclicvoltammetry'

def test_bagit_converter():
    target = target_dir + source_dir + '/File053_BagIt.zip'
    with tempfile.TemporaryDirectory() as td:
      with zipfile.ZipFile(target, 'r') as z:
            z.extractall(td)
        
      bgcv = BagItBaseConverter(td)
    assert len(bgcv.data) == 3
    assert len(bgcv.images) == 3