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

def test_base_composer_generate_original_metadata_headers(jcamp_file_1h):
    base_converter = JcampBaseConverter(jcamp_file_1h)
    composer = BaseComposer(core=base_converter)
    headers = composer._BaseComposer__header_original_metadata()
    assert headers == [
      '\n',
      '$$ === CHEMSPECTRA ORIGINAL METADATA ===\n',
    ]

def test_base_composer_generate_original_metadata(jcamp_file_1h):
    base_converter = JcampBaseConverter(jcamp_file_1h)
    metadata = {"__comments": ["first comment, second comment"], ".AVERAGES": "16", ".SHIFT REFERENCE": ["INTERNAL", "CDCl3", 1, 15.91938], "#.AVERAGES": "16",}
    base_converter.dic = metadata
    composer = BaseComposer(core=base_converter)
    auto_metadata = composer.generate_original_metadata()
    assert auto_metadata == [
      '\n',
      '$$ === CHEMSPECTRA ORIGINAL METADATA ===\n',
      '###.AVERAGES= 16\n',
      '###.SHIFT REFERENCE= INTERNAL, CDCl3, 1, 15.91938\n',
      '###.AVERAGES= 16\n',
      '\n\n'
    ]

def test_base_composer_generate_original_metadata_keeps_xypoints_when_not_ntuples(jcamp_file_1h):
    # XYPOINTS is the authoritative data table for XYPOINTS-classified
    # spectra (e.g. HPLC) and must still appear in the original-metadata
    # dump when the source isn't NTUPLES.
    base_converter = JcampBaseConverter(jcamp_file_1h)
    metadata = {
        "DATACLASS": ["XYPOINTS"],
        ".AVERAGES": "16",
        "XYPOINTS": ["(XY..XY)\n999, 999;\n888, 888;"],
    }
    base_converter.dic = metadata
    composer = BaseComposer(core=base_converter)
    auto_metadata = composer.generate_original_metadata()
    assert auto_metadata == [
      '\n',
      '$$ === CHEMSPECTRA ORIGINAL METADATA ===\n',
      '###.AVERAGES= 16\n',
      '###XYPOINTS= (XY..XY)\n999, 999;\n888, 888;\n',
      '\n\n'
    ]

def test_base_composer_generate_original_metadata_ignores_xypoints_for_ntuples(jcamp_file_1h):
    # A stray top-level XYPOINTS LDR alongside NTUPLES-classified data (e.g.
    # MASS SPECTRUM / LC-MS) duplicates the real NTUPLES data table and must
    # not leak into the CHEMSPECTRA ORIGINAL METADATA dump, same as
    # XYDATA/XYDATA_OLD already don't.
    base_converter = JcampBaseConverter(jcamp_file_1h)
    metadata = {
        "DATACLASS": ["NTUPLES"],
        ".AVERAGES": "16",
        "XYPOINTS": ["(XY..XY)\n999, 999;\n888, 888;"],
        "XYDATA_OLD": ["(XY..XY)\n999, 999;\n888, 888;"],
    }
    base_converter.dic = metadata
    composer = BaseComposer(core=base_converter)
    auto_metadata = composer.generate_original_metadata()
    assert auto_metadata == [
      '\n',
      '$$ === CHEMSPECTRA ORIGINAL METADATA ===\n',
      '###.AVERAGES= 16\n',
      '\n\n'
    ]