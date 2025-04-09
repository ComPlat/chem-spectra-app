from werkzeug.datastructures import FileStorage

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter
from chem_spectra.lib.composer.ni import NIComposer
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
    nicv = JcampNIConverter(jbcv)
    assert nicv.is_cyclic_volta == True
    assert nicv.datatype == 'CYCLIC VOLTAMMETRY'
    assert nicv.xs[0] == 1.49048
    assert nicv.ys[0] == 5.34724E-06


def test_cv_compose():
    target = target_dir + source_dir + '/RCV_LSH-R444_full+Fc.jdx'
    jbcv = JcampBaseConverter(target)
    nicv = JcampNIConverter(jbcv)
    nicmpsr = NIComposer(nicv)
    assert "$$ === CHEMSPECTRA CYCLIC VOLTAMMETRY ===\n" in nicmpsr.meta
    assert "##$CSSCANRATE=0.09\n" in nicmpsr.meta
    assert "##$CSSPECTRUMDIRECTION=NEGATIVE\n" in nicmpsr.meta
    assert "##$CSCYCLICVOLTAMMETRYDATA=\n" in nicmpsr.meta
