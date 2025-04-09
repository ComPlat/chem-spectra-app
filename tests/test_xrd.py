from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter

target_dir = './tests/fixtures/'
source_dir = 'source/'


def test_xrd_converter_1():
    target = target_dir + source_dir + '/xrd/Test_data_op_001_XRD.jdx'
    jbcv = JcampBaseConverter(target)
    assert jbcv.is_xrd == True


def test_xrd_converter_2():
    target = target_dir + source_dir + '/xrd/Test_data_op_002_XRD.jdx'
    jbcv = JcampBaseConverter(target)
    assert jbcv.is_xrd == True


def test_is_non_nmr():
    file1 = target_dir + source_dir + '/xrd/Test_data_op_001_XRD.jdx'
    file2 = target_dir + source_dir + '/xrd/Test_data_op_002_XRD.jdx'
    jbcv1 = JcampBaseConverter(file1)
    jbcv2 = JcampBaseConverter(file2)
    assert jbcv1.non_nmr == True
    assert jbcv2.non_nmr == True
