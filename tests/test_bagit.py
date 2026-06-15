import tempfile
import zipfile

from chem_spectra.lib.converter.bagit.base import BagItBaseConverter

from tests.dataset_catalog import dataset_path_str

def test_bagit_converter():
    target = dataset_path_str('CV-B-001')
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(target, 'r') as z:
            z.extractall(td)
          
        bgcv = BagItBaseConverter(td)
    assert len(bgcv.data) == 3
    assert len(bgcv.images) == 3
