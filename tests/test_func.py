import json
import io
from werkzeug.datastructures import FileStorage

from chem_spectra.model.molecule import MoleculeModel

target_dir = './tests/fixtures/'
source_dir = 'source/'


def test_molfile2chem():
    with open(target_dir + source_dir + '/molfile/svs813f1_B.mol', 'rb') as f:
        molfile = FileStorage(f).stream.read().decode('utf-8')
        mm = MoleculeModel(molfile)

    assert mm.can() == 'CCC1CCC(=O)C(C21SCCS2)C'
    assert mm.mass() == 230.079907196
