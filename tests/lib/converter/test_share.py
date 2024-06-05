import numpy as np
import pytest
from chem_spectra.lib.converter.share import (
  parse_params, parse_solvent, reduce_pts
)

@pytest.fixture
def expected_default_params():
    default_itg = {'stack': [], 'refArea': 1, 'refFactor': 1, 'shift': 0}
    default_mpy = {'stack': [], 'smExtext': False, 'shift': 0}
    default_wavelength = {'name': 'CuKalpha', 'value': 0.15406, 'label': 'Cu K-alpha', 'unit': 'nm'}
    return {
        'select_x': None,
        'ref_name': None,
        'ref_value': None,
        'peaks_str': None,
        'delta': 0.0,
        'mass': 0,
        'scan': None,
        'thres': None,
        'clear': False,
        'integration': default_itg,
        'multiplicity': default_mpy,
        'fname': '',
        'waveLength': default_wavelength,
        'list_max_min_peaks': None,
        'cyclicvolta': None,
        'jcamp_idx': 0,
        'axesUnits': None,
        'detector': None,
        'dsc_meta_data': None,
    }

def test_parse_params_without_params(expected_default_params):
    parsed_data = parse_params(None)
    assert expected_default_params == parsed_data
    
def test_parse_params_select_x():
    params = {'select_x': None}
    parsed_data = parse_params(params)
    assert parsed_data['select_x'] is None
    
    params = {'select_x': 0.0}
    parsed_data = parse_params(params)
    assert parsed_data['select_x'] == 0.0
    
def test_parse_params_ref_name():
    params = {'ref_name': None}
    parsed_data = parse_params(params)
    assert parsed_data['ref_name'] is None
    
    params = {'ref_name': 'just a text'}
    parsed_data = parse_params(params)
    assert parsed_data['ref_name'] == 'just a text'

def test_parse_params_ref_value():
    params = {'ref_value': None}
    parsed_data = parse_params(params)
    assert parsed_data['ref_value'] is None
    
    params = {'ref_value': 0.0}
    parsed_data = parse_params(params)
    assert parsed_data['ref_value'] == 0.0
    
def test_parse_params_peaks_str():
    params = {'peaks_str': None}
    parsed_data = parse_params(params)
    assert parsed_data['peaks_str'] is None
    
    params = {'peaks_str': 'just a text'}
    parsed_data = parse_params(params)
    assert parsed_data['peaks_str'] == 'just a text'

def test_parse_params_delta():
    params = {'select_x': None}
    parsed_data = parse_params(params)
    assert parsed_data['delta'] == 0.0
    
    params = {'select_x': 0.0}
    parsed_data = parse_params(params)
    assert parsed_data['delta'] == 0.0
    
    params = {'select_x': 1.0, 'ref_name': '- - -'}
    parsed_data = parse_params(params)
    assert parsed_data['delta'] == 0.0
    
    params = {'select_x': 1.0, 'ref_name': 'just a text'}
    parsed_data = parse_params(params)
    assert parsed_data['delta'] == 0.0
    
    params = {'select_x': 1.0, 'ref_name': 'just a text', 'ref_value': 1.5}
    parsed_data = parse_params(params)
    assert parsed_data['delta'] == 0.5
    
def test_parse_params_mass():
    params = {'mass': None}
    parsed_data = parse_params(params)
    assert parsed_data['mass'] == 0
    
    params = {'mass': 1.5}
    parsed_data = parse_params(params)
    assert parsed_data['mass'] == 1.5

def test_parse_params_scan():
    params = {'scan': None}
    parsed_data = parse_params(params)
    assert parsed_data['scan'] is None
    
    params = {'scan': 2}
    parsed_data = parse_params(params)
    assert parsed_data['scan'] == 2

def test_parse_params_thres():
    params = {'thres': None}
    parsed_data = parse_params(params)
    assert parsed_data['thres'] is None
    
    params = {'thres': 2}
    parsed_data = parse_params(params)
    assert parsed_data['thres'] == 2

def test_parse_params_clear():
    params = {'clear': None}
    parsed_data = parse_params(params)
    assert parsed_data['clear'] == False
    
    params = {'clear': False}
    parsed_data = parse_params(params)
    assert parsed_data['clear'] == False
    
    params = {'clear': True}
    parsed_data = parse_params(params)
    assert parsed_data['clear'] == True

def test_parse_params_ext():
    params = {'ext': None}
    parsed_data = parse_params(params)
    assert parsed_data['ext'] == ''
    
    params = {'ext': '.jdx'}
    parsed_data = parse_params(params)
    assert parsed_data['ext'] == '.jdx'

def test_parse_params_fname():
    params = {'fname': 'just a normal text'}
    parsed_data = parse_params(params)
    assert parsed_data['fname'] == ''
    
    params = {'fname': 'original.jdx'}
    parsed_data = parse_params(params)
    assert parsed_data['fname'] == 'original'
    
    params = {'fname': 'original.changed.jdx'}
    parsed_data = parse_params(params)
    assert parsed_data['fname'] == 'original.changed'
    
    params = {'fname': 'original.peak.jdx'}
    parsed_data = parse_params(params)
    assert parsed_data['fname'] == 'original'
    
    params = {'fname': 'original.edit.jdx'}
    parsed_data = parse_params(params)
    assert parsed_data['fname'] == 'original'

def test_parse_params_integration():
    #TODO: need to be updated
    assert 1==1

def test_parse_params_multiplicity():
    #TODO: need to be updated
    assert 1==1

def test_parse_params_waveLength():
    #TODO: need to be updated
    assert 1==1

def test_parse_params_list_max_min_peaks():
    #TODO: need to be updated
    assert 1==1
    
def test_parse_params_cyclic_volta():
    #TODO: need to be updated
    assert 1==1
    
def test_parse_params_jcamp_idx():
    params = {'jcamp_idx': None}
    parsed_data = parse_params(params)
    assert parsed_data['jcamp_idx'] == 0
    
    params = {'jcamp_idx': 1}
    parsed_data = parse_params(params)
    assert parsed_data['jcamp_idx'] == 1
    
def test_parse_params_axesUnits():
    params = {'axesUnits': None}
    parsed_data = parse_params(params)
    assert parsed_data['axesUnits'] is None
    
    params = {'axesUnits': '{"axes": null}'}
    parsed_data = parse_params(params)
    assert parsed_data['axesUnits'] is None
    
    params = {'axesUnits': '{"axes": []}'}
    parsed_data = parse_params(params)
    assert parsed_data['axesUnits'] is None
    
    params = {'jcamp_idx': 1, 'axesUnits': '{"axes": [{"xUnit": "label x", "yUnit": "label y"}]}'}
    parsed_data = parse_params(params)
    assert parsed_data['axesUnits'] is None

    params = {'axesUnits': '{"axes": [{"xUnit": "label x", "yUnit": "label y"}]}'}
    parsed_data = parse_params(params)
    assert parsed_data['axesUnits'] == {"xUnit": "label x", "yUnit": "label y"}

def test_parse_params_data_type_mapping():
    #TODO: need to be updated
    assert 1==1
    
def test_parse_params_detector():
    #TODO: need to be updated
    assert 1==1

def test_parse_solvent():
    #TODO: need to be updated
    assert 1==1

def test_parse_dsc_meta_data():
    params = {'dsc_meta_data': '{"meltingPoint": "1.0", "tg": "1.0"}'}
    parsed_data = parse_params(params)
    assert parsed_data['dsc_meta_data'] == {"meltingPoint": "1.0", "tg": "1.0"}
    
def test_reduce_pts_when_does_not_have_any_x():
    array_data = []
    data = np.array(array_data)
    reduced_data = reduce_pts(data)
    assert np.array_equal(reduced_data, data)

def test_reduce_pts_when_under_limitation():
    array_data = [[x, x+1] for x in range(0, 3999)]
    data = np.array(array_data)
    reduced_data = reduce_pts(data)
    assert np.array_equal(reduced_data, data)
    
def test_reduce_pts_when_reached_limitation():
    array_data = [[x, x+1] for x in range(0, 4500)]
    data = np.array(array_data)
    expected_array_data = [[x, x+1] for x in range(576, 4500)]
    expected_data = np.array(expected_array_data)
    reduced_data = reduce_pts(data)
    assert np.array_equal(reduced_data, expected_data)
    