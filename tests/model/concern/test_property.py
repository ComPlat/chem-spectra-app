import pytest
from werkzeug.datastructures import FileStorage
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.model.concern.property import decorate_sim_property, __simulate_nmr

target_dir = './tests/fixtures/'
source_dir = 'source/'

class TestConverterObject:
    def __init__(self, layout='1H'):
        self.ncl = layout
        self.simu_peaks = False

@pytest.fixture
def parameters():
    normal_1H_jbcv = TestConverterObject()
    normal_13C_jbcv = TestConverterObject(layout='13C')


    return {"normal_1H_jbcv": normal_1H_jbcv, "normal_13C_jbcv": normal_13C_jbcv}

def test_simulate_nmr_1H_without_molfile(parameters):
    normal_1H_jbcv = parameters["normal_1H_jbcv"]
    jbcv = __simulate_nmr(normal_1H_jbcv, None)
    assert jbcv["invalid_molfile"] is True
    assert jbcv["origin_jbcv"] is not None
    assert jbcv["origin_jbcv"]  == normal_1H_jbcv

def test_simulate_nmr_13C_without_molfile(parameters):
    normal_13C_jbcv = parameters["normal_13C_jbcv"]
    jbcv = __simulate_nmr(normal_13C_jbcv, None)
    assert jbcv["invalid_molfile"] is True
    assert jbcv["origin_jbcv"] is not None
    assert jbcv["origin_jbcv"]  == normal_13C_jbcv

def test_simulate_nmr_with_invalid_molfile(parameters):
    with open(target_dir + source_dir + '/molfile/invalid_molfile.mol', 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    normal_1H_jbcv = parameters["normal_1H_jbcv"]
    jbcv = __simulate_nmr(normal_1H_jbcv, molfile)
    assert jbcv["invalid_molfile"] is True
    assert jbcv["origin_jbcv"] is not None
    assert jbcv["origin_jbcv"]  == normal_1H_jbcv

def test_simulate_nmr_1H_valid_molfile(parameters):
    with open(target_dir + source_dir + '/molfile/svs813f1_B.mol', 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    normal_1H_jbcv = parameters["normal_1H_jbcv"]

    jbcv = __simulate_nmr(normal_1H_jbcv, molfile)
    assert jbcv is not None
    assert isinstance(jbcv, TestConverterObject)

def test_simulate_nmr_1H_valid_molfile(parameters):
    with open(target_dir + source_dir + '/molfile/svs813f1_B.mol', 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    normal_13C_jbcv = parameters["normal_13C_jbcv"]

    jbcv = __simulate_nmr(normal_13C_jbcv, molfile)
    assert jbcv is not None
    assert isinstance(jbcv, TestConverterObject)

def test_decorate_sim_property_no_nmrsimulate_no_molfile(parameters):
    normal_1H_jbcv = parameters["normal_1H_jbcv"]

    jbcv = decorate_sim_property(normal_1H_jbcv, None)
    assert jbcv == normal_1H_jbcv

def test_decorate_sim_property_1H_no_nmrsimulate(parameters):
    with open(target_dir + source_dir + '/molfile/svs813f1_B.mol', 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    normal_1H_jbcv = parameters["normal_1H_jbcv"]

    jbcv = decorate_sim_property(normal_1H_jbcv, molfile)
    assert jbcv == normal_1H_jbcv


def test_decorate_sim_property_13C_no_nmrsimulate(parameters):
    with open(target_dir + source_dir + '/molfile/svs813f1_B.mol', 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    normal_13C_jbcv = parameters["normal_13C_jbcv"]

    jbcv = decorate_sim_property(normal_13C_jbcv, molfile)
    assert jbcv == normal_13C_jbcv

def test_decorate_sim_property_no_molfile(parameters):
    normal_1H_jbcv = parameters["normal_1H_jbcv"]

    jbcv = decorate_sim_property(normal_1H_jbcv, None, True)
    assert jbcv == normal_1H_jbcv

def test_decorate_sim_property_1H_invalid_molfile(parameters):
    with open(target_dir + source_dir + '/molfile/invalid_molfile.mol', 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    normal_1H_jbcv = parameters["normal_1H_jbcv"]

    jbcv = decorate_sim_property(normal_1H_jbcv, molfile, True)
    assert jbcv is not None
    assert jbcv["invalid_molfile"] is True
    assert jbcv["origin_jbcv"] is not None
    assert jbcv["origin_jbcv"]  == normal_1H_jbcv

def test_decorate_sim_property_1H(parameters):
    with open(target_dir + source_dir + '/molfile/svs813f1_B.mol', 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    normal_1H_jbcv = parameters["normal_1H_jbcv"]

    jbcv = decorate_sim_property(normal_1H_jbcv, molfile, True)
    assert jbcv is not None

def test_decorate_sim_property_13C_invalid_molfile(parameters):
    with open(target_dir + source_dir + '/molfile/invalid_molfile.mol', 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    normal_13C_jbcv = parameters["normal_13C_jbcv"]

    jbcv = decorate_sim_property(normal_13C_jbcv, molfile, True)
    assert jbcv is not None
    assert jbcv["invalid_molfile"] is True
    assert jbcv["origin_jbcv"] is not None
    assert jbcv["origin_jbcv"]  == normal_13C_jbcv

def test_decorate_sim_property_13C(parameters):
    with open(target_dir + source_dir + '/molfile/svs813f1_B.mol', 'rb') as f:
        molfile = FileContainer(FileStorage(f))

    normal_13C_jbcv = parameters["normal_13C_jbcv"]

    jbcv = decorate_sim_property(normal_13C_jbcv, molfile, True)
    assert jbcv is not None
