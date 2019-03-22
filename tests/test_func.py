import json
import io
from werkzeug.datastructures import FileStorage

from fixtures.mock_predict import RequestPredictNmr
from chem_spectra.model.converter.chem import molfile2chem
from chem_spectra.model.converter.ms import MsConverter
from chem_spectra.model.composer.ms import MsComposer

target_dir = './tests/fixtures/'
source_dir = 'source/'


def test_molfile2chem():
    with open(target_dir + source_dir + '/molfile/svs813f1_B.mol', 'rb') as f:
        molfile = FileStorage(f)
        smi, mass = molfile2chem(molfile)

    assert smi == 'CCC1CCC(=O)C(C21SCCS2)C'
    assert mass == 230.079907196


def test_ms_converter_composer():
    params = { 'mass': 230.079907196 }

    with open(target_dir + source_dir + '/mzml/svs813f1.mzML', 'rb') as f:
        file = FileStorage(f)
        mcv = MsConverter(file, params)
        mcp = MsComposer(mcv)

    lines = mcp.tf_jcamp().read()[:800].decode('utf-8').split('\n')

    assert '##$SCANAUTOTARGET=3' in lines
    assert '##$SCANEDITTARGET=' in lines
    assert '##$SCANCOUNT=24' in lines
    assert '51.012176513671875, 34359.0' in lines
