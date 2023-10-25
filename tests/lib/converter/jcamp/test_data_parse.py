import pytest
from chem_spectra.lib.converter.jcamp.data_parse import read_parsed_jdx_data

@pytest.fixture
def parsed_data_with_old_dic_structure():
    return {'DATATYPE': ['NMR']}, []

@pytest.fixture
def parsed_data_with_new_dic_structure():
    return {
      '_datatype_':[
        {'DATATYPE': ['NMR']}
      ]
    }, []

@pytest.fixture
def parsed_data_with_ntupes_type_xy_xy():
    return {
      '_datatype_':[
        {
          'DATATYPE': ['NMR'],
          'DATACLASS': ['NTUPLES'],
          'DATATABLE': ['(XY..XY), PEAKS\n51.012176513671875, 34359.0']
        }
      ]
    }, []

def test_read_parsed_jdx_data_with_old_dic(parsed_data_with_old_dic_structure):
    dic, data = read_parsed_jdx_data(parsed_data_with_old_dic_structure)
    assert isinstance(dic, dict)
    assert dic['DATATYPE'] == ['NMR']
    assert isinstance(data, list)

def test_read_parsed_jdx_data_with_new_dic(parsed_data_with_new_dic_structure):
    dic, data = read_parsed_jdx_data(parsed_data_with_new_dic_structure)
    assert isinstance(dic, dict)
    assert dic['DATATYPE'] == ['NMR']
    assert isinstance(data, list)

def test_read_parsed_jdx_data_with_ntupes_xy_xy(parsed_data_with_ntupes_type_xy_xy):
    dic, data = read_parsed_jdx_data(parsed_data_with_ntupes_type_xy_xy)
    assert isinstance(dic, dict)
    assert dic['DATATYPE'] == ['NMR']
    assert isinstance(data, list)
