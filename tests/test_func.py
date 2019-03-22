import json
import io
from werkzeug.datastructures import FileStorage
from chem_spectra.model.converter.chem import molfile2chem
from fixtures.mock_predict import RequestPredictNmr


target_dir = './tests/fixtures/'
source_dir = 'source/'


def test_molfile2chem():
    with open(target_dir + source_dir + '/molfile/svs813f1_B.mol', 'rb') as f:
        molfile = FileStorage(f)
        smi, mass = molfile2chem(molfile)

    assert smi == 'CCC1CCC(=O)C(C21SCCS2)C'
    assert mass == 230.079907196
