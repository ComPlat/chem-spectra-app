import numpy as np
import pytest
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter
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
        'lcms_uvvis_wavelength': None,
        'lcms_mz_page': None,
        'lcms_mz_page_data': None,
        'transmittance': False,
        'invert_y': None,
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


def test_parse_solvent_uses_existing_ref_name_for_1h_solvent_ranges():
    base = type('MockBase', (), {
        'ncl': '1H',
        'params': {'ref_name': None},
        'dic': {
            '$CSSOLVENTNAME': ['DMSO-d6 (quin)'],
            '$CSSOLVENTVALUE': ['2.5'],
            '$CSSOLVENTX': ['0'],
            '.SOLVENTNAME': [''],
            '.SHIFTREFERENCE': [''],
        },
        'solv_peaks': [],
    })()

    parse_solvent(base)

    assert base.solv_peaks == [(2.48, 2.52)]
    assert base.dic['$CSSOLVENTNAME'] == ['DMSO-d6 (quin)']


def test_parse_solvent_uses_existing_ref_name_for_13c_solvent_ranges():
    base = type('MockBase', (), {
        'ncl': '13C',
        'params': {'ref_name': None},
        'dic': {
            '$CSSOLVENTNAME': ['Chloroform-d (t)'],
            '$CSSOLVENTVALUE': ['77.16'],
            '$CSSOLVENTX': ['0'],
            '.SOLVENTNAME': [''],
            '.SHIFTREFERENCE': [''],
        },
        'solv_peaks': [],
    })()

    parse_solvent(base)

    assert base.solv_peaks == [(74.16, 80.16)]
    assert base.dic['$CSSOLVENTNAME'] == ['Chloroform-d (t)']


def test_parse_solvent_falls_back_when_existing_ref_is_auto_offset():
    base = type('MockBase', (), {
        'ncl': '1H',
        'params': {'ref_name': None},
        'dic': {
            '$CSSOLVENTNAME': ['AUTO-OFFSET'],
            '$CSSOLVENTVALUE': ['0.000000'],
            '$CSSOLVENTX': ['4.238162'],
            '.SOLVENTNAME': ['DMSO'],
            '.SHIFTREFERENCE': ['DMSO'],
        },
        'solv_peaks': [],
    })()

    parse_solvent(base)

    assert base.solv_peaks == [(2.48, 2.52)]
    assert base.dic['$CSSOLVENTNAME'] == ['DMSO-d6 (quin)']
    assert base.dic['$CSSOLVENTVALUE'] == ['2.5']
    assert base.dic['$CSSOLVENTX'] == ['0']


def test_parse_solvent_an_explicit_unsupported_reference_turns_filtering_off():
    # the request's reference outranks the file's, even when it names no
    # solvent with known ranges
    base = type('MockBase', (), {
        'ncl': '1H',
        'params': {'ref_name': 'TMS'},
        'dic': {
            '$CSSOLVENTNAME': ['DMSO-d6 (quin)'],
            '$CSSOLVENTVALUE': ['2.5'],
            '$CSSOLVENTX': ['0'],
        },
        'solv_peaks': [],
    })()

    parse_solvent(base)

    assert base.solv_peaks == []
    assert base.dic['$CSSOLVENTNAME'] == ['DMSO-d6 (quin)']


def test_parse_solvent_auto_offset_without_a_bruker_solvent():
    # an acqus without SOLVENT leaves both fields None
    base = type('MockBase', (), {
        'ncl': '1H',
        'params': {'ref_name': None},
        'dic': {
            '$CSSOLVENTNAME': ['AUTO-OFFSET'],
            '$CSSOLVENTVALUE': ['0.000000'],
            '$CSSOLVENTX': ['4.009459'],
            '.SOLVENTNAME': None,
            '.SHIFTREFERENCE': None,
        },
        'solv_peaks': [],
    })()

    parse_solvent(base)

    assert base.solv_peaks == []
    assert base.dic['$CSSOLVENTNAME'] == ['AUTO-OFFSET']


def test_parse_dsc_meta_data():
    params = {'dsc_meta_data': '{"meltingPoint": "1.0", "tg": "1.0"}'}
    parsed_data = parse_params(params)
    assert parsed_data['dsc_meta_data'] == {"meltingPoint": "1.0", "tg": "1.0"}


def test_parse_lcms_mz_page_data():
    params = {'lcms_mz_page_data': '{"0":[{"x":101.1,"y":1000}]}'}
    parsed_data = parse_params(params)
    assert parsed_data['lcms_mz_page_data'] == '{"0":[{"x":101.1,"y":1000}]}'
    
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
    

def test_parse_params_tolerates_a_partial_cyclicvolta_payload():
    """The ELN does not guarantee spectraList or the entry's 'list'.

    ViewSpectra.js reads `spectraList?.[curveIdx]` and bails when it is
    missing, so a payload without it is not a client error. Subscripting it
    unconditionally made this a 500.
    """
    import json
    from chem_spectra.lib.converter.share import parse_params

    assert parse_params(
        {'cyclic_volta': json.dumps({'cvDisplay': {}})}
    )['list_max_min_peaks'] is None

    assert parse_params(
        {'cyclic_volta': json.dumps({'spectraList': [{}]})}
    )['list_max_min_peaks'] is None

    assert parse_params({
        'cyclic_volta': json.dumps({'spectraList': [{'list': []}]}),
        'jcamp_idx': 5,
    })['list_max_min_peaks'] is None


def test_cv_shift_offset_tolerates_a_spectra_entry_without_shift():
    from chem_spectra.lib.shared.calc import (
        cal_cyclic_volta_shift_prev_offset_at_index,
    )

    assert cal_cyclic_volta_shift_prev_offset_at_index(
        {'spectraList': [{}]}, 0) == 0.0
    assert cal_cyclic_volta_shift_prev_offset_at_index(
        {'spectraList': [{'hasRefPeak': True, 'shift': {'prevValue': 2.0}}]},
        0) == 2.0


@pytest.mark.parametrize('value, expected', [
    (None, None), ('', None), ('  ', None),
    (True, True), ('true', True), ('TRUE', True), ('1', True), ('yes', True),
    (1, True),
    (False, False), ('false', False), (' False ', False), ('0', False),
    ('no', False), (0, False),
    # unrecognised: the same as not sent, never an explicit false
    ('undefined', None), ('null', None), ('on', None), ('t', None),
    ('y', None), ('off', None), (2, None),
])
def test_invert_y_is_tristate(value, expected):
    """An explicit false clears ##$CSINVERTY, so only a recognised false may
    produce one. `undefined` is what JS FormData sends for a missing value."""
    assert parse_params({'invert_y': value})['invert_y'] is expected


# The solvents react-spectra-editor offers as a reference (`list_shift.js`)
# that this table knows, with the editor's value for each signal. A name the
# table does not match turns solvent filtering off, and a signal outside every
# window is left among the sample peaks.
EDITOR_SOLVENT_SIGNALS = [
    ('13C', 'Acetone-d6 (sep)', 29.640),
    ('13C', 'Acetone-d6 (broad)', 206.260),
    ('13C', 'Acetonitrile-d3 (sep)', 1.32),
    ('13C', 'Acetonitrile-d3 (s)', 118.26),
    ('13C', 'Benzene (t)', 128.060),
    ('13C', 'Chloroform-d (t)', 77.16),
    ('13C', 'Dichloromethane-d2 (quin)', 53.84),
    ('13C', 'DMSO-d6', 39.52),
    ('13C', 'Methanol-d4 (sep)', 49.00),
    ('1H', 'Acetone-d6 (quin)', 2.05),
    ('1H', 'Acetonitrile-d3 (quin)', 1.94),
    ('1H', 'Benzene (s)', 7.16),
    ('1H', 'Chloroform-d (s)', 7.26),
    ('1H', 'D2O (s)', 4.79),
    ('1H', 'Dichloromethane-d2 (t)', 5.32),
    ('1H', 'DMSO-d6 (quin)', 2.50),
    ('1H', 'Methanol-d4 (quin)', 3.31),
    ('1H', 'Methanol-d4 (s)', 4.87),
]


@pytest.mark.parametrize('ncl, name, value', EDITOR_SOLVENT_SIGNALS)
def test_every_signal_the_editor_offers_has_a_window(ncl, name, value):
    from chem_spectra.lib.converter.share import _solvent_config
    config = _solvent_config(ncl, name)
    assert config is not None, name
    _, _, windows = config
    assert any(low < value < high for low, high in windows), (name, windows)


def test_the_acetonitrile_nitrile_line_is_collapsed_too():
    # CD3 septet at 1.32 and the nitrile carbon at 118.26: one line kept for
    # each, and a sample peak outside both windows kept as well
    converter = JcampTechniqueConverter.__new__(JcampTechniqueConverter)
    converter.ncl = '13C'
    converter.obs_freq = 100.6
    converter.solv_peaks = [(-1.68, 4.32), (115.26, 121.26)]
    step = 21.0 / 100.6
    septet = [{'x': 1.32 + k * step, 'y': w}
              for k, w in zip(range(-3, 4), [1.0, 3.0, 6.0, 7.0, 6.0, 3.0, 1.0])]
    nitrile = [{'x': 118.26, 'y': 4.0}]
    sample = [{'x': 125.4, 'y': 9.0}]

    filtered = converter._JcampTechniqueConverter__filter_solvent_peaks(
        septet + nitrile + sample)

    assert sorted(round(p['x'], 2) for p in filtered) == [1.32, 118.26, 125.4]
