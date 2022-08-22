from werkzeug.datastructures import FileStorage

from chem_spectra.controller.helper.file_container import FileContainer

from chem_spectra.lib.converter.nmrium.base import NMRiumDataConverter

source = './tests/fixtures/source/nmrium/nmrium_test.nmrium'

def test_file_is_none():
    nmrium_converter = NMRiumDataConverter(None)
    assert nmrium_converter.data == None

def test_read_file_content():
    with open(source, 'rb') as f:
        nmriumFile = FileContainer(FileStorage(f))

    nmrium_converter = NMRiumDataConverter(nmriumFile)
    assert nmrium_converter.data != None
    assert len(nmrium_converter.data['x']) == 65536