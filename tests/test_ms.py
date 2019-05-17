from werkzeug.datastructures import FileStorage

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ms import JcampMSConverter
from chem_spectra.lib.converter.ms import MSConverter
from chem_spectra.lib.composer.ms import MSComposer
from chem_spectra.controller.helper.file_container import FileContainer

target_dir = './tests/fixtures/'
source_dir = 'source/'


def test_ms_mzml_converter_composer():
    params = {'mass': 230.079907196}

    with open(target_dir + source_dir + '/ms/svs813f1.mzML', 'rb') as f:
        file = FileStorage(f)
        file = FileContainer(file)
        mscv = MSConverter(file, params)
        mscp = MSComposer(mscv)

    lines = mscp.tf_jcamp().read()[:800] \
                .decode('utf-8', errors='ignore').split('\n')

    assert '##$SCANAUTOTARGET=3' in lines
    assert '##$SCANEDITTARGET=3' in lines
    assert '##$SCANCOUNT=24' in lines
    assert '##$THRESHOLD=0.05' in lines
    assert '51.012176513671875, 34359.0' in lines


def test_ms_jcamp_converter_composer():
    params = {'mass': 230.079907196}

    target = target_dir + source_dir + '/ms/svs813f1.jdx'
    jbcv = JcampBaseConverter(target, params)
    mscv = JcampMSConverter(jbcv)
    mscp = MSComposer(mscv)

    lines = mscp.tf_jcamp().read()[:800] \
                .decode('utf-8', errors='ignore').split('\n')

    assert '##$SCANAUTOTARGET=3' in lines
    assert '##$SCANEDITTARGET=16' in lines
    assert '##$SCANCOUNT=24' in lines
    assert '##$THRESHOLD=0.155' in lines
    assert '51.012176513671875, 34359.0' in lines
