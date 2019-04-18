import json
import io
from werkzeug.datastructures import FileStorage

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.model.molecule import MoleculeModel

target_dir = './tests/fixtures/'
source_dir = 'source/'


def test_molfile2chem():
    with open(target_dir + source_dir + '/molfile/svs813f1_B.mol', 'rb') as f:
        molfile = FileContainer(FileStorage(f))
        mm = MoleculeModel(molfile)

    assert mm.smi == 'CCC1CCC(=O)C(C)C12SCCS2'
    assert mm.mass == 230.079907196
