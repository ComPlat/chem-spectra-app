from werkzeug.datastructures import FileStorage

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.model.molecule import MoleculeModel

from tests.dataset_catalog import dataset_path


def test_molfile2chem():
    with dataset_path('MOL-002').open('rb') as f:
        molfile = FileContainer(FileStorage(f))
        mm = MoleculeModel(molfile)

    assert mm.smi == 'CCC1CCC(=O)C(C)C12SCCS2'
    assert mm.mass == 230.079907196
