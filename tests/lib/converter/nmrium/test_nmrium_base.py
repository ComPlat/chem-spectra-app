import pytest
from werkzeug.datastructures import FileStorage
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.lib.converter.nmrium.base import NMRiumDataConverter

source_files = [
    './tests/fixtures/source/nmrium/nmrium_test_v3.nmrium',
    './tests/fixtures/source/nmrium/nmrium_test_v4.nmrium',
    './tests/fixtures/source/nmrium/nmrium_test_v5.nmrium'
]


@pytest.fixture
def nmrium_files():
    files = []
    for source in source_files:
        with open(source, 'rb') as f:
            nmriumFile = FileContainer(FileStorage(f))
            files.append(nmriumFile)
    return files


@pytest.fixture
def parsed_json_nmrium_data():
    dic_data = {
        "spectra": [
            {
                'id': 'id1',
                'data': {'x': [1, 2, 3], 're': [4, 5, 6]},
                'info': {'isFid': True}
            },
            {
                'id': 'id2',
                'data': {'x': [1.1, 2.1, 3.1], 're': [4.1, 5.1, 6.1]},
                'info': {'isFid': False}
            }
        ],
        "correlations": {
            "values": [
                {
                    'link': [
                        {
                            'experimentID': 'id2'
                        }
                    ]
                }
            ]
        }
    }
    return dic_data


def test_file_is_none():
    nmrium_converter = NMRiumDataConverter(None)
    assert nmrium_converter.data == None


def test_find_displaying_spectrum_failed():
    nmrium_converter = NMRiumDataConverter()
    spectrum_data = nmrium_converter._NMRiumDataConverter__find_displaying_spectra()
    assert spectrum_data == None


def test_find_displaying_spectrum(parsed_json_nmrium_data):
    nmrium_converter = NMRiumDataConverter()
    spectrum_data = nmrium_converter._NMRiumDataConverter__find_displaying_spectra(
        parsed_json_nmrium_data['spectra'])
    assert spectrum_data != None
    assert spectrum_data == [parsed_json_nmrium_data['spectra'][1]]


def test_find_displaying_spectrum_id_failed():
    nmrium_converter = NMRiumDataConverter()
    spectrum_id = nmrium_converter._NMRiumDataConverter__find_displaying_spectrum_id()
    assert spectrum_id == ''


def test_find_displaying_spectrum_id(parsed_json_nmrium_data):
    nmrium_converter = NMRiumDataConverter()
    spectrum_id = nmrium_converter._NMRiumDataConverter__find_displaying_spectrum_id(
        parsed_json_nmrium_data['correlations'])
    assert spectrum_id == 'id2'


def test_parsing_spectra_none():
    nmrium_converter = NMRiumDataConverter()
    spectrum_data = nmrium_converter._NMRiumDataConverter__parsing_spectra()
    assert spectrum_data == None


def test_parsing_spectra(parsed_json_nmrium_data):
    expected_value = parsed_json_nmrium_data['spectra'][1]
    nmrium_converter = NMRiumDataConverter()
    spectrum_data = nmrium_converter._NMRiumDataConverter__parsing_spectra(
        parsed_json_nmrium_data)
    assert spectrum_data != None
    assert spectrum_data == expected_value


def test_read_file(nmrium_files):
    for file in nmrium_files:
        nmrium_converter = NMRiumDataConverter()
        nmrium_converter.file = file
        spectrum_data = nmrium_converter._NMRiumDataConverter__read_file()
        assert spectrum_data is not None


def test_read_file_content(nmrium_files):
    for file in nmrium_files:
        nmrium_converter = NMRiumDataConverter(file)
        assert nmrium_converter.data != None
        assert len(nmrium_converter.data['x']) == 65536
        assert len(nmrium_converter.data['x']) == len(
            nmrium_converter.data['y'])


def test_read_xy_values_failed():
    nmrium_converter = NMRiumDataConverter()
    x_values, y_values = nmrium_converter._NMRiumDataConverter__read_xy_values(
        None)
    assert x_values == None
    assert y_values == None


def test_read_xy_values(parsed_json_nmrium_data):
    expected_x_values = [1.1, 2.1, 3.1]
    expected_y_values = [4.1, 5.1, 6.1]

    nmrium_converter = NMRiumDataConverter()
    spectrum_data = nmrium_converter._NMRiumDataConverter__parsing_spectra(
        parsed_json_nmrium_data)
    x_values, y_values = nmrium_converter._NMRiumDataConverter__read_xy_values(
        spectrum_data)
    assert x_values != None
    assert y_values != None
    assert len(x_values) == len(x_values)
    assert x_values == expected_x_values
    assert y_values == expected_y_values


def test_parsing_xy_values_failed():
    nmrium_converter = NMRiumDataConverter()
    xy_values = nmrium_converter._NMRiumDataConverter__parsing_xy_values(None)
    assert xy_values == None


def test_parsing_xy_values(parsed_json_nmrium_data):
    expected_x_values = [1.1, 2.1, 3.1]
    expected_y_values = [4.1, 5.1, 6.1]

    nmrium_converter = NMRiumDataConverter()
    spectrum_data = nmrium_converter._NMRiumDataConverter__parsing_spectra(
        parsed_json_nmrium_data)
    xy_values = nmrium_converter._NMRiumDataConverter__parsing_xy_values(
        spectrum_data)
    assert xy_values != None
    assert xy_values['x'] == expected_x_values
    assert xy_values['y'] == expected_y_values
