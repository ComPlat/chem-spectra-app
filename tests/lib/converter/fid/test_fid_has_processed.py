import tempfile
import zipfile
import os

from chem_spectra.lib.converter.fid.base import FidBaseConverter
from chem_spectra.lib.converter.fid.bruker import FidHasBruckerProcessed

source_dir_1h = './tests/fixtures/source/bruker/1H.zip'
filename_1h = '1H'
params_1h = {'fname':'1H.zip'}

source_dir_13c = './tests/fixtures/source/bruker/13C.zip'
filename_13c = '13C'
params_13c = {'fname':'13C.zip'}

def test_read_data():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(source_dir_1h, 'r') as z:
            z.extractall(td)
        
        target_dir = os.path.join(td, '1')
        fid_conveter = FidHasBruckerProcessed(target_dir=target_dir, params=params_1h, fname=filename_1h)
        assert len(fid_conveter.data) == 2

def test_has_fid_data():
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(source_dir_1h, 'r') as z:
            z.extractall(td)
        
        target_dir = os.path.join(td, '1')
        fid_conveter = FidHasBruckerProcessed(target_dir=target_dir, params=params_1h, fname=filename_1h)
        assert type(fid_conveter.data[0]) is FidBaseConverter
    