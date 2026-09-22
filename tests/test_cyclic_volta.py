from werkzeug.datastructures import FileStorage

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter
from chem_spectra.lib.composer.technique import TechniqueComposer
from chem_spectra.controller.helper.file_container import FileContainer

target_dir = './tests/fixtures/'
source_dir = 'source/cyclicvoltammetry'

def test_cv_base_converter():
    target = target_dir + source_dir + '/RCV_LSH-R444_full+Fc.jdx'
    jbcv = JcampBaseConverter(target)
    assert jbcv.is_cyclic_volta == True

def test_cv_ni_converter():
    target = target_dir + source_dir + '/RCV_LSH-R444_full+Fc.jdx'
    jbcv = JcampBaseConverter(target)
    tcv = JcampTechniqueConverter(jbcv)
    assert tcv.technique.cyclic_voltammetry is True
    assert tcv.datatype == 'CYCLIC VOLTAMMETRY'
    assert tcv.xs[0] == 1.49048
    assert tcv.ys[0] == 5.34724E-06

def test_cv_compose():
    target = target_dir + source_dir + '/RCV_LSH-R444_full+Fc.jdx'
    jbcv = JcampBaseConverter(target)
    tcv = JcampTechniqueConverter(jbcv)
    nicmpsr = TechniqueComposer(tcv)
    assert "$$ === CHEMSPECTRA CYCLIC VOLTAMMETRY ===\n" in nicmpsr.meta
    assert "##$CSSCANRATE=0.09\n" in nicmpsr.meta
    assert "##$CSSPECTRUMDIRECTION=NEGATIVE\n" in nicmpsr.meta
    assert "##$CSCYCLICVOLTAMMETRYDATA=\n" in nicmpsr.meta
