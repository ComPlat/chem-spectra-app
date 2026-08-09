from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter
from chem_spectra.lib.composer.ni import NIComposer

from tests.dataset_catalog import dataset_path_str

CV_FIXTURE = dataset_path_str('CV-009')

def test_cv_base_converter():
    jbcv = JcampBaseConverter(CV_FIXTURE)
    assert jbcv.is_cyclic_volta == True

def test_cv_ni_converter():
    jbcv = JcampBaseConverter(CV_FIXTURE)
    nicv = JcampNIConverter(jbcv)
    assert nicv.is_cyclic_volta == True
    assert nicv.datatype == 'CYCLIC VOLTAMMETRY'
    assert nicv.xs[0] == 1.49048
    assert nicv.ys[0] == 5.34724E-06

def test_cv_compose():
    jbcv = JcampBaseConverter(CV_FIXTURE)
    nicv = JcampNIConverter(jbcv)
    nicmpsr = NIComposer(nicv)
    assert "$$ === CHEMSPECTRA CYCLIC VOLTAMMETRY ===\n" in nicmpsr.meta
    assert "##$CSSCANRATE=0.09\n" in nicmpsr.meta
    assert "##$CSSPECTRUMDIRECTION=NEGATIVE\n" in nicmpsr.meta
    assert "##$CSCYCLICVOLTAMMETRYDATA=\n" in nicmpsr.meta
