"""Classification of a JCAMP file onto a spectrum kind.

Covers the three defects fixed together: an unrecognised ``##DATA TYPE=``
used to raise ``UnboundLocalError`` (HTTP 500), block selection preferred
whatever ``data_type.json`` happened to list last, and classification
preferred whatever it happened to list first.
"""

import json

import pytest

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter

source_nmr = './tests/fixtures/source/1H.dx'
source_ir = './tests/fixtures/source/IR.dx'
source_hplc = './tests/fixtures/source/hplc/chromatogram.jdx'
source_chi = './tests/fixtures/source/CHI-224_10.jdx'

HPLC_HEADER = '##DATA TYPE=HPLC UV/VIS SPECTRUM'


@pytest.fixture
def jcamp_with_datatype(tmp_path):
    """Build a parsable JCAMP whose primary ##DATA TYPE= is `datatype`."""
    template = open(source_hplc).read()

    def _build(datatype):
        target = tmp_path / 'probe.jdx'
        target.write_text(template.replace(HPLC_HEADER, '##DATA TYPE=' + datatype, 1))
        return str(target)

    return _build


# - - - unrecognised datatypes - - -

UNMAPPED = ['SQUID', 'TENSIOMETRY', 'LINEAR SWEEP VOLTAMMETRY',
            'SINGLE CRYSTAL X-RAY DIFFRACTION', 'NEUTRON SCATTERING']


@pytest.mark.parametrize('datatype', UNMAPPED)
def test_unrecognised_datatype_does_not_raise(jcamp_with_datatype, datatype):
    # these raised UnboundLocalError from __index_target and surfaced as a 500
    converter = JcampTechniqueConverter(JcampBaseConverter(jcamp_with_datatype(datatype)))
    # nothing is recognised, so the file's only block is read as the curve
    assert converter.target.datatype == datatype


@pytest.mark.parametrize('datatype', UNMAPPED)
def test_unrecognised_datatype_is_not_treated_as_nmr(jcamp_with_datatype, datatype):
    # non_nmr False would give it chemical-shift axes and multiplet analysis
    base = JcampBaseConverter(jcamp_with_datatype(datatype))
    assert base.typ == ''
    assert base.technique.key != 'NMR'


def test_unrecognised_datatype_is_logged(jcamp_with_datatype, caplog):
    JcampBaseConverter(jcamp_with_datatype('SQUID'), {'fname': 'probe.jdx'})
    assert 'SQUID' in caplog.text
    # parse_params strips the extension, so the log carries the stem
    assert 'probe' in caplog.text


# - - - newly mapped aliases - - -

@pytest.mark.parametrize('datatype,expected', [
    ('CONTINUOUS MASS SPECTRUM', 'MS'),
    ('GEL PERMEATION CHROMATOGRAPHY', 'SIZE EXCLUSION CHROMATOGRAPHY'),
])
def test_alias_maps_onto_existing_kind(jcamp_with_datatype, datatype, expected):
    assert JcampBaseConverter(jcamp_with_datatype(datatype)).typ == expected


# - - - ordering - - -

def test_auxiliary_blocks_stay_unmapped():
    """NMR FID and peak tables must not be in data_type.json.

    Block selection takes the first *recognised* datatype, so mapping an
    auxiliary block would make every Bruker-derived file read that block
    instead of the spectrum.
    """
    import json
    import os
    from chem_spectra.lib.converter.jcamp import base as base_module

    path = os.path.join(os.path.dirname(base_module.__file__), 'data_type.json')
    with open(path) as handle:
        mapped = {v.upper() for vals in json.load(handle)['datatypes'].values() for v in vals}

    # the first five occur in this repo's own fixtures (PEAK ASSIGNMENTS in 12
    # blocks, NMR FID in 8, NMRPEAKTABLE and NMR PEAK ASSIGNMENTS in 2 each,
    # NMR PEAK TABLE in 1); the rest are spellings chemotion-converter-app can
    # emit, including its malformed 'NMP' variant
    for auxiliary in ['PEAK ASSIGNMENTS', 'NMR FID', 'NMRPEAKTABLE',
                      'NMR PEAK ASSIGNMENTS', 'NMR PEAK TABLE',
                      'NMP PEAK ASSIGNMENTS', 'INFRARED PEAK TABLE',
                      'INFRARED INTERFEROGRAM']:
        assert auxiliary not in mapped


def test_file_order_decides_between_two_recognised_datatypes(tmp_path):
    """MASS TIC precedes MASS SPECTRUM in the file, so the file is LC/MS.

    Previously the key order of data_type.json decided, and MS (key 3)
    outranked LC/MS (key 17) regardless of what the file said.
    """
    template = open(source_hplc).read()
    target = tmp_path / 'tic_then_ms.jdx'
    target.write_text(
        template.replace(HPLC_HEADER, '##DATA TYPE=MASS TIC\n##DATA TYPE=MASS SPECTRUM', 1)
    )
    assert JcampBaseConverter(str(target)).typ == 'LC/MS'


@pytest.mark.parametrize('path,expected_typ,expected_datatype', [
    (source_nmr, 'NMR', 'NMR SPECTRUM'),
    (source_ir, 'INFRARED', 'INFRARED SPECTRUM'),
    (source_hplc, 'HPLC UVVIS', 'HPLC UV/VIS SPECTRUM'),
])
def test_known_files_classify_and_select_block_unchanged(
        path, expected_typ, expected_datatype):
    """Named by what the chosen block declares rather than by its position.

    The old assertion was an index into a merged LDR list, which is the thing
    the block-aware read removes: for `1H.dx` that index was 1 while the block
    is the file's third.
    """
    base = JcampBaseConverter(path)
    assert base.typ == expected_typ
    assert JcampTechniqueConverter(base).target.datatype == expected_datatype


# - - - caller-supplied mapping - - -
#
# Callers send this as `data_type_mapping`; parse_params renames it to
# `user_data_type_mapping` internally. Passing the internal name raises
# nothing -- the built-in mapping is silently used instead -- so these tests
# go through the caller-facing name on purpose.

USER_MAPPING = json.dumps({'datatypes': {
    'NMR': ['NMR SPECTRUM', 'NMRSPECTRUM'],
    'SQUID': ['SQUID', 'SQUID MAGNETOMETRY'],
    'X-RAY DIFFRACTION': ['X-RAY DIFFRACTION'],
}})


@pytest.fixture
def user_mapping_params():
    return {'data_type_mapping': USER_MAPPING}


@pytest.mark.parametrize('datatype', ['SQUID', 'SQUID MAGNETOMETRY'])
def test_user_mapping_recognises_its_own_datatypes(
        jcamp_with_datatype, user_mapping_params, datatype):
    # SQUID is absent from the built-in data_type.json; the caller adds it,
    # and both the key and its alias must resolve to the key
    base = JcampBaseConverter(jcamp_with_datatype(datatype), user_mapping_params)
    assert base.typ == 'SQUID'
    assert base.technique.key != 'NMR'


def test_user_mapping_replaces_the_builtin_one(user_mapping_params):
    # HPLC UVVIS is in the built-in mapping but not in the caller's, so it
    # must become unrecognised rather than falling back. This combination
    # used to raise UnboundLocalError from __index_target.
    base = JcampBaseConverter(source_hplc, user_mapping_params)
    assert base.typ == ''
    assert base.technique.key != 'NMR'
    # unrecognised under this mapping, so the first block carrying data is read
    assert JcampTechniqueConverter(base).target.datatype == 'HPLC UV/VIS SPECTRUM'


def test_user_mapping_still_classifies_nmr(user_mapping_params):
    base = JcampBaseConverter(source_nmr, user_mapping_params)
    converter = JcampTechniqueConverter(base)
    assert base.typ == 'NMR'
    assert base.technique.key == 'NMR'
    assert converter.target.datatype == 'NMR SPECTRUM'
    assert converter.threshold == 0.005


def test_user_mapping_drives_the_threshold(jcamp_with_datatype, user_mapping_params):
    # __thres() keeps its own copy of the mapping load, so it is checked
    # separately from classification
    base = JcampBaseConverter(jcamp_with_datatype('X-RAY DIFFRACTION'), user_mapping_params)
    assert base.typ == 'X-RAY DIFFRACTION'
    assert JcampTechniqueConverter(base).threshold == 1.00


@pytest.mark.parametrize('value', ['', None])
def test_absent_user_mapping_falls_back_to_builtin(value):
    base = JcampBaseConverter(source_hplc, {'data_type_mapping': value})
    assert base.typ == 'HPLC UVVIS'
    assert JcampTechniqueConverter(base).threshold == 0.05


def test_example_mapping_stays_in_sync_with_the_live_one():
    """`data_type.json.example` seeds a fresh install and had drifted.

    `spectra_layout_api.load_data_types()` copies it into place on
    FileNotFoundError, so a new deploy comes up with whatever it contains.
    It was missing the entire LC/MS key before this branch; without this
    test the same drift recurs while the suite stays green.
    """
    import json
    import os
    from chem_spectra.lib.converter.jcamp import base as base_module

    directory = os.path.dirname(base_module.__file__)
    with open(os.path.join(directory, 'data_type.json')) as handle:
        live = json.load(handle)
    with open(os.path.join(directory, 'data_type.json.example')) as handle:
        example = json.load(handle)

    assert example == live


def test_block_selection_agrees_with_classification(tmp_path):
    """Classification and block selection must not point at different blocks.

    CHI-224_10.jdx is LINK / NMR SPECTRUM / NMR PEAK TABLE. Renaming the
    trailing auxiliary block to MASS TIC makes it a datatype that is mapped
    but sits last in the flattened value order of data_type.json, while
    NMR SPECTRUM sits first. Before #291 the selection took the last
    match and chose the MASS TIC block while
    __set_datatype classified the file as NMR from the first matching key --
    the two disagreed. Both now take the first recognised block in the
    file's own order.
    """
    body = open(source_chi).read().replace(
        '##DATA TYPE=\tNMR PEAK TABLE', '##DATA TYPE=\tMASS TIC', 1,
    )
    target = tmp_path / 'link_nmr_then_tic.jdx'
    target.write_text(body)

    base = JcampBaseConverter(str(target))
    assert base.datatypes == ['LINK', 'NMR SPECTRUM', 'MASS TIC']
    assert base.typ == 'NMR'
    # the NMR SPECTRUM block, not the MASS TIC one that classification did
    # not choose -- now stated as the block rather than as an index
    assert JcampTechniqueConverter(base).target.datatype == 'NMR SPECTRUM'


def test_warning_points_at_the_mapping_that_is_actually_in_effect(
        jcamp_with_datatype, user_mapping_params, caplog):
    """A caller-supplied mapping replaces the built-in one.

    Telling such a caller to 'add it to data_type.json' is useless advice --
    that file is not consulted for their request.
    """
    JcampBaseConverter(jcamp_with_datatype('NEUTRON SCATTERING'), user_mapping_params)
    assert 'data_type_mapping supplied with this request' in caplog.text
    assert 'data_type.json' not in caplog.text

    caplog.clear()
    JcampBaseConverter(jcamp_with_datatype('NEUTRON SCATTERING'))
    assert 'data_type.json' in caplog.text
