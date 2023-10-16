import pytest
from chem_spectra.lib.composer.base import BaseComposer
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter

source_nmr = './tests/fixtures/source/1H.dx'

@pytest.fixture
def jcamp_file_1h():
    return source_nmr  

def test_init_ni_composer_success(jcamp_file_1h):
    base_converter = JcampBaseConverter(jcamp_file_1h)
    composer = BaseComposer(core=base_converter)
    
    assert composer is not None
    assert composer.core == base_converter

def test_base_composer_generate_auto_metadata_headers(jcamp_file_1h):
    base_converter = JcampBaseConverter(jcamp_file_1h)
    composer = BaseComposer(core=base_converter)
    headers = composer._BaseComposer__header_auto_metadata()
    assert headers == [
      '\n',
      '$$ === CHEMSPECTRA AUTO METADATA ===\n',
      '##$CSAUTOMETADATA=\n',
    ]

def test_base_composer_generate_auto_metadata(jcamp_file_1h):
    base_converter = JcampBaseConverter(jcamp_file_1h)
    metadata = {"just a string": "just a string value"}
    base_converter.auto_metadata = metadata
    composer = BaseComposer(core=base_converter)
    auto_metadata = composer.generate_auto_metadata()
    assert auto_metadata == [
      '\n',
      '$$ === CHEMSPECTRA AUTO METADATA ===\n',
      '##$CSAUTOMETADATA=\n',
      'JUST A STRING=just a string value\n'
    ]
