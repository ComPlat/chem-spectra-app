from werkzeug.datastructures import FileStorage

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ms import JcampMSConverter
from chem_spectra.lib.converter.ms import MSConverter
from chem_spectra.lib.composer.ms import MSComposer
from chem_spectra.controller.helper.file_container import FileContainer

from tests.dataset_catalog import dataset_path_str

def test_xrd_converter_1():
    target = dataset_path_str('XRD-004')
    jbcv = JcampBaseConverter(target)
    assert jbcv.is_xrd == True

def test_xrd_converter_2():
    target = dataset_path_str('XRD-005')
    jbcv = JcampBaseConverter(target)
    assert jbcv.is_xrd == True

def test_is_non_nmr():
    file1 = dataset_path_str('XRD-004')
    file2 = dataset_path_str('XRD-005')
    jbcv1 = JcampBaseConverter(file1)
    jbcv2 = JcampBaseConverter(file2)
    assert jbcv1.non_nmr == True
    assert jbcv2.non_nmr == True
