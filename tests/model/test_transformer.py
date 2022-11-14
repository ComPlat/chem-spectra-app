import tempfile
import zipfile
import os

from chem_spectra.lib.converter.fid.base import FidBaseConverter
from chem_spectra.lib.converter.fid.bruker import FidHasBruckerProcessed
from chem_spectra.model.transformer import TransformerModel

source_dir_1h = './tests/fixtures/source/bruker/1H.zip'
filename_1h = '1H'
params_1h = {'fname':'1H.zip'}

source_dir_molfile = './tests/fixtures/source/molfile/c60h57fn4.mol'

def test_zip2cv_with_processed_file():
    file = open(source_dir_1h, "rb")
    molfile = open(source_dir_molfile, "r")
    tranform_model = TransformerModel(file, molfile=molfile, params=params_1h)
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(source_dir_1h, 'r') as z:
            z.extractall(td)
        
        target_dir = os.path.join(td, '1')

        list_converters, list_composers = tranform_model.zip2cv_with_processed_file(target_dir=target_dir, params=params_1h, file_name=filename_1h)
        assert len(list_converters) == 2
        assert len(list_composers) == 2
        