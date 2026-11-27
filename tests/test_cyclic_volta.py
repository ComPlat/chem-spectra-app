from werkzeug.datastructures import FileStorage

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.curve import JcampCurveConverter
from chem_spectra.lib.composer.curve import CurveComposer
from chem_spectra.controller.helper.file_container import FileContainer

target_dir = './tests/fixtures/'
source_dir = 'source/cyclicvoltammetry'

def test_cv_base_converter():
    target = target_dir + source_dir + '/RCV_LSH-R444_full+Fc.jdx'
    jbcv = JcampBaseConverter(target)
    assert jbcv.technique.cyclic_voltammetry is True

def test_cv_ni_converter():
    target = target_dir + source_dir + '/RCV_LSH-R444_full+Fc.jdx'
    jbcv = JcampBaseConverter(target)
    ccv = JcampCurveConverter(jbcv)
    assert ccv.technique.cyclic_voltammetry is True
    assert ccv.datatype == 'CYCLIC VOLTAMMETRY'
    assert ccv.xs[0] == 1.49048
    assert ccv.ys[0] == 5.34724E-06

def test_cv_compose():
    target = target_dir + source_dir + '/RCV_LSH-R444_full+Fc.jdx'
    jbcv = JcampBaseConverter(target)
    ccv = JcampCurveConverter(jbcv)
    nicmpsr = CurveComposer(ccv)
    assert "$$ === CHEMSPECTRA CYCLIC VOLTAMMETRY ===\n" in nicmpsr.meta
    assert "##$CSSCANRATE=0.09\n" in nicmpsr.meta
    assert "##$CSSPECTRUMDIRECTION=NEGATIVE\n" in nicmpsr.meta
    assert "##$CSCYCLICVOLTAMMETRYDATA=\n" in nicmpsr.meta
