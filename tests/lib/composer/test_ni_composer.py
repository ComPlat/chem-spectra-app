import json
import pytest
import numpy as np
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter
from chem_spectra.lib.composer.ni import NIComposer
from chem_spectra.controller.helper.file_container import FileContainer

source_nmr = './tests/fixtures/source/1H.dx'
source_nmr_edit = './tests/fixtures/source/1H.edit.jdx'
source_ir = './tests/fixtures/source/IR.dx'
source_dir_molfile = './tests/fixtures/source/molfile/c60h57fn4.mol'

@pytest.fixture
def data_schema():
  return [
        './tests/fixtures/source/nmrium/nmrium_data_schema_v3.json',
        './tests/fixtures/source/nmrium/nmrium_data_schema_v4.json',
        './tests/fixtures/source/nmrium/nmrium_data_schema_v5.json'
    ]

@pytest.fixture
def jcamp_file_1h():
    return source_nmr

@pytest.fixture
def jcamp_file_1h_edit():
    return source_nmr_edit

@pytest.fixture
def jcamp_file_ir():
    return source_ir

@pytest.fixture
def molfile_data():
    molfile = open(source_dir_molfile, "r")
    molfile_str = molfile.read()
    return molfile_str

def test_init_ni_composer_failed():
    with pytest.raises(Exception) as error:
        _ = NIComposer(None)
        
    assert error is not None    

def test_init_ni_composer_success(jcamp_file_1h):
    base_converter = JcampBaseConverter(jcamp_file_1h)
    ni_converter = JcampNIConverter(base=base_converter)
    ni_composer = NIComposer(core=ni_converter)
    
    assert ni_composer is not None
    assert ni_composer.core == ni_converter

def test_ni_composer_original_metadata(jcamp_file_1h):
    base_converter = JcampBaseConverter(jcamp_file_1h)
    ni_converter = JcampNIConverter(base=base_converter)
    ni_composer = NIComposer(core=ni_converter)

    assert ni_composer is not None
    assert '$$ === CHEMSPECTRA ORIGINAL METADATA ===\n' in ni_composer.meta

def test_ni_composer_header(jcamp_file_1h):
    base_converter = JcampBaseConverter(jcamp_file_1h)
    ni_converter = JcampNIConverter(base=base_converter)
    ni_composer = NIComposer(core=ni_converter)
    headers = ni_composer._NIComposer__header_base()
    assert headers == [
        '\n',
        '$$ === CHEMSPECTRA SPECTRUM ORIG ===\n',
        '##TITLE=\n',
        '##JCAMP-DX=5.00\n',
        '##DATA TYPE=NMR SPECTRUM\n',
        '##DATA CLASS=XYDATA\n',
        '##$CSCATEGORY=SPECTRUM\n',
        '##ORIGIN=PH\n',
        '##OWNER=PH\n'
    ]

def test_ni_composer_generate_nmrium_not_nmr(jcamp_file_ir):
    base_converter = JcampBaseConverter(jcamp_file_ir)
    ni_converter = JcampNIConverter(base=base_converter)
    ni_composer = NIComposer(core=ni_converter)
    nmrium_file = ni_composer.generate_nmrium()
    
    assert nmrium_file is None

def test_ni_composer_check_schema_v3(data_schema, jcamp_file_1h):
    schema_source = data_schema[0]
    schema = ''
    with open(schema_source, 'r') as f:
        schema = f.read()
    dic_schema = json.loads(schema)
    expected_schema_keys = np.array(dic_schema.keys())
    
    base_converter = JcampBaseConverter(jcamp_file_1h)
    ni_converter = JcampNIConverter(base=base_converter)
    ni_composer = NIComposer(core=ni_converter)
    nmrium_file = ni_composer.generate_nmrium(version=3)
    with open(nmrium_file.name, 'r') as f:
        nmrium = f.read()
        nmrium_spectra = json.loads(nmrium)
        nmrium_spectra_keys = np.array(nmrium_spectra.keys())
        assert np.array_equal(expected_schema_keys, nmrium_spectra_keys)
        assert nmrium_spectra['version'] == 3

def test_ni_composer_check_schema_v4(data_schema, jcamp_file_1h):
    schema_source = data_schema[1]
    schema = ''
    with open(schema_source, 'r') as f:
        schema = f.read()
    dic_schema = json.loads(schema)
    expected_schema_keys = np.array(dic_schema.keys())
    
    base_converter = JcampBaseConverter(jcamp_file_1h)
    ni_converter = JcampNIConverter(base=base_converter)
    ni_composer = NIComposer(core=ni_converter)
    nmrium_file = ni_composer.generate_nmrium()
    with open(nmrium_file.name, 'r') as f:
        nmrium = f.read()
        nmrium_spectra = json.loads(nmrium)
        nmrium_spectra_keys = np.array(nmrium_spectra.keys())
        assert np.array_equal(expected_schema_keys, nmrium_spectra_keys)
        assert nmrium_spectra['version'] == 4

def test_ni_composer_check_schema_v5(data_schema, jcamp_file_1h):
    schema_source = data_schema[2]
    schema = ''
    with open(schema_source, 'r') as f:
        schema = f.read()
    dic_schema = json.loads(schema)
    expected_schema_keys = np.array(dic_schema.keys())
    
    base_converter = JcampBaseConverter(jcamp_file_1h)
    ni_converter = JcampNIConverter(base=base_converter)
    ni_composer = NIComposer(core=ni_converter)
    nmrium_file = ni_composer.generate_nmrium(version=5)
    with open(nmrium_file.name, 'r') as f:
        nmrium = f.read()
        nmrium_spectra = json.loads(nmrium)
        nmrium_spectra_keys = np.array(nmrium_spectra.keys())
        assert np.array_equal(expected_schema_keys, nmrium_spectra_keys)
        assert nmrium_spectra['version'] == 5

def test_ni_composer_nmrium_spectra(jcamp_file_1h):
    base_converter = JcampBaseConverter(jcamp_file_1h)
    ni_converter = JcampNIConverter(base=base_converter)
    ni_composer = NIComposer(core=ni_converter)
    nmrium_spectra = ni_composer._NIComposer__generate_nmrim_spectra()
    
    assert len(nmrium_spectra) == 1

    xy_points = nmrium_spectra[0]['data']
    x_values = xy_points['x']
    y_values = xy_points['re']
    assert len(x_values) == 65536
    assert len(y_values) == 65536

def test_ni_composer_nmrium_peaks(jcamp_file_1h):
    base_converter = JcampBaseConverter(jcamp_file_1h)
    ni_converter = JcampNIConverter(base=base_converter)
    ni_composer = NIComposer(core=ni_converter)
    nmrium_peaks = ni_composer._NIComposer__generate_nmrim_peaks()
    
    assert 'values' in nmrium_peaks
    assert 'options' in nmrium_peaks

    peaks = nmrium_peaks['values']
    assert len(peaks) == 46

def test_ni_composer_nmrium_integrals(jcamp_file_1h_edit):
    base_converter = JcampBaseConverter(jcamp_file_1h_edit)
    ni_converter = JcampNIConverter(base=base_converter)
    ni_composer = NIComposer(core=ni_converter)
    nmrium_integrals = ni_composer._NIComposer__generate_nmrim_integrals()
    
    assert 'values' in nmrium_integrals
    assert 'options' in nmrium_integrals

    integral = nmrium_integrals['values']
    assert len(integral) == 1

def test_ni_composer_nmrium_ranges(jcamp_file_1h_edit):
    base_converter = JcampBaseConverter(jcamp_file_1h_edit)
    ni_converter = JcampNIConverter(base=base_converter)
    ni_composer = NIComposer(core=ni_converter)
    nmrium_integrals = ni_composer._NIComposer__generate_nmrim_ranges()
    
    assert 'values' in nmrium_integrals
    assert 'options' in nmrium_integrals

    integral = nmrium_integrals['values']
    assert len(integral) == 1

def test_ni_composer_generate_nmrium(jcamp_file_1h):
    base_converter = JcampBaseConverter(jcamp_file_1h)
    ni_converter = JcampNIConverter(base=base_converter)
    ni_composer = NIComposer(core=ni_converter)
    nmrium_file = ni_composer.generate_nmrium()
    
    assert nmrium_file is not None
    with open(nmrium_file.name) as file:
        file_content = file.read()
        assert file_content != ""

        dic_content = json.loads(file_content)
        dic_data = dic_content['data']
        assert dic_data['spectra'] is not None

        xy_points = dic_data['spectra'][0]['data']
        x_values = xy_points['x']
        y_value = xy_points['re']
        assert len(x_values) == 65536
        assert len(y_value) == 65536

def test_ni_composer_generate_nmrium_without_molfile(jcamp_file_1h):
    base_converter = JcampBaseConverter(jcamp_file_1h)
    ni_converter = JcampNIConverter(base=base_converter)
    ni_composer = NIComposer(core=ni_converter)
    nmrium_file = ni_composer.generate_nmrium()
    
    assert nmrium_file is not None
    with open(nmrium_file.name) as file:
        file_content = file.read()
        assert file_content != ""

        dic_content = json.loads(file_content)
        dic_data = dic_content['data']
        molecules = dic_data['molecules']
        assert molecules is not None
        assert len(molecules) == 0

def test_ni_composer_generate_nmrium_with_molfile(jcamp_file_1h, molfile_data):
    base_converter = JcampBaseConverter(jcamp_file_1h)
    ni_converter = JcampNIConverter(base=base_converter)
    ni_composer = NIComposer(core=ni_converter)
    nmrium_file = ni_composer.generate_nmrium(molfile_data=molfile_data)
    
    assert nmrium_file is not None
    with open(nmrium_file.name) as file:
        file_content = file.read()
        assert file_content != ""

        dic_content = json.loads(file_content)
        dic_data = dic_content['data']
        molecules = dic_data['molecules']
        assert molecules is not None
        assert len(molecules) == 1

        molecule_data = molecules[0]
        assert 'id' in molecule_data.keys()
        assert molecule_data['label'] == 'P1'
        assert molecule_data['molfile'] == molfile_data
