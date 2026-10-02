"""An unrecognised ##DATA TYPE= must survive into the composed file.

Since #291 a file whose `##DATA TYPE=` is not in `data_type.json` no longer
crashes -- it takes the generic curve path. But `__set_datatype` returned `''`
for it, the composer writes that value straight back out
(`composer/technique.py`), and `DATATYPE` is on the suppression list for the
original-metadata dump (`composer/base.py`). Between them, the only record of
what the file said it was got erased at compose time.

The spectrum rendered fine either way, so nothing failed. What was lost is the
ability to reclassify the file later -- which is precisely what happens when an
under-specified technique (SQUID, TENSIOMETRY, LSV ...) is added to
`data_type.json` after the fact: every file already processed has forgotten
what it was.

`SQUID` is used here because it is a real `##DATA TYPE=` emitted by
chemotion-converter-app that this app deliberately does not map yet.
"""

import re

import pytest

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter
from chem_spectra.lib.composer.technique import TechniqueComposer

RECOGNISED = './tests/fixtures/source/IR.dx'


def write_with_datatype(tmp_path, datatype, name='probe.dx'):
    """A real fixture with only its ##DATA TYPE= swapped.

    Building the file from a real one rather than by hand keeps every other
    LDR consistent, so a failure here means the datatype handling changed and
    not that the fixture stopped parsing.
    """
    source = open(RECOGNISED).read()
    swapped = re.sub(
        r'##DATA TYPE=.*INFRARED SPECTRUM', f'##DATA TYPE={datatype}',
        source, flags=re.IGNORECASE,
    )
    assert swapped != source, 'the fixture no longer has the header to swap'
    target = tmp_path / name
    target.write_text(swapped)
    return str(target)


def compose(path):
    core = JcampBaseConverter(path)
    composer = TechniqueComposer(JcampTechniqueConverter(core))
    lines = composer.tf_jcamp()
    body = lines.read() if hasattr(lines, 'read') else ''.join(lines)
    return body.decode() if isinstance(body, bytes) else body


def datatypes_in(text):
    return [m.strip() for m in re.findall(r'##DATA TYPE=(.*)', text)]


def test_an_unrecognised_datatype_is_kept_on_the_converter(tmp_path):
    core = JcampBaseConverter(write_with_datatype(tmp_path, 'SQUID'))
    assert core.datatype == 'SQUID'


def test_it_is_still_not_a_known_technique(tmp_path):
    """Preserving the label must not smuggle it past classification."""
    core = JcampBaseConverter(write_with_datatype(tmp_path, 'SQUID'))
    assert core.typ == ''


def test_it_reaches_the_composed_file(tmp_path):
    composed = compose(write_with_datatype(tmp_path, 'SQUID'))
    assert 'SQUID' in datatypes_in(composed)


def test_it_survives_a_round_trip(tmp_path):
    """The case that matters: the stored file can still be reclassified."""
    first = tmp_path / 'first.jdx'
    first.write_text(compose(write_with_datatype(tmp_path, 'SQUID')))
    assert JcampBaseConverter(str(first)).datatype == 'SQUID'
    second = compose(str(first))
    assert 'SQUID' in datatypes_in(second)


def test_structural_blocks_are_never_mistaken_for_the_measurement(tmp_path):
    """A composed file whose spectrum block lost its datatype (anything
    written before this fix) holds only LINK and PEAKTABLE markers. None of
    them names a technique, so the answer stays empty rather than 'LINK'."""
    core = JcampBaseConverter(write_with_datatype(tmp_path, 'SQUID'))
    assert core._JcampBaseConverter__unrecognised_datatype() == 'SQUID'
    core.datatypes = ['LINK', 'PEAKTABLE', 'PEAKTABLE']
    assert core._JcampBaseConverter__unrecognised_datatype() == ''
    core.datatypes = ['LINK', 'NMR FID']
    assert core._JcampBaseConverter__unrecognised_datatype() == ''


# The spellings `test_auxiliary_blocks_stay_unmapped` pins as non-measurements.
# Kept in step with that test deliberately: it guarantees none of these is in
# data_type.json, which is what sends a file carrying one down this fallback in
# the first place.
AUXILIARY_SPELLINGS = [
    'PEAK ASSIGNMENTS', 'NMR FID', 'NMRPEAKTABLE', 'NMR PEAK ASSIGNMENTS',
    'NMR PEAK TABLE', 'NMP PEAK ASSIGNMENTS', 'INFRARED PEAK TABLE',
    'INFRARED INTERFEROGRAM',
]


@pytest.mark.parametrize('auxiliary', AUXILIARY_SPELLINGS)
def test_an_auxiliary_block_never_wins_over_the_measurement(tmp_path, auxiliary):
    """The fallback takes the *first* non-structural datatype, so a derived
    block sitting ahead of the real one would be preserved in its place.

    The first version of this only skipped `LINK`, `NMR FID` and anything
    ending in `PEAKTABLE` -- which is how this app composes the block, but not
    how chemotion-converter-app emits it (`NMR PEAK TABLE`, with spaces). Six
    of these eight spellings slipped through and were returned as though they
    named the measurement.
    """
    core = JcampBaseConverter(write_with_datatype(tmp_path, 'SQUID'))
    core.datatypes = ['LINK', auxiliary, 'SQUID']
    assert core._JcampBaseConverter__unrecognised_datatype() == 'SQUID'


@pytest.mark.parametrize('auxiliary', AUXILIARY_SPELLINGS)
def test_every_auxiliary_spelling_is_recognised_as_one(auxiliary):
    assert JcampBaseConverter._is_auxiliary_datatype(auxiliary)


@pytest.mark.parametrize('measurement', [
    'SQUID', 'TENSIOMETRY', 'LINEAR SWEEP VOLTAMMETRY',
    'SINGLE CRYSTAL X-RAY DIFFRACTION', 'INFRARED TRANSFERED SPECTRUM',
])
def test_no_real_measurement_is_mistaken_for_auxiliary(measurement):
    """The exclusion must not grow teeth. Every unmapped entry in the
    converter app's dropdown that names a real measurement stays preserved."""
    assert not JcampBaseConverter._is_auxiliary_datatype(measurement)


@pytest.mark.parametrize('datatype,typ', [
    ('X-RAY DIFFRACTION', 'X-RAY DIFFRACTION'),
    ('RAMAN SPECTRUM', 'RAMAN'),
    ('CYCLIC VOLTAMMETRY', 'CYCLIC VOLTAMMETRY'),
])
def test_recognised_files_still_classify(tmp_path, datatype, typ):
    """The fallback must only ever run when nothing matched."""
    core = JcampBaseConverter(write_with_datatype(tmp_path, datatype))
    assert core.datatype == datatype
    assert core.typ == typ


def test_the_untouched_fixture_is_unaffected():
    """The real IR fixture, with no swap at all."""
    core = JcampBaseConverter(RECOGNISED)
    assert core.datatype == 'INFRARED SPECTRUM'
    assert core.typ == 'INFRARED'


def test_single_crystal_xrd_is_the_one_cross_repo_hazard(tmp_path):
    """Preserving the real string removes an accidental shield. Read this
    before shipping.

    react-spectra-editor's `readLayout` matches `##DATA TYPE=` by *substring*;
    this backend matches exactly. `SINGLE CRYSTAL X-RAY DIFFRACTION` contains
    `X-RAY DIFFRACTION`, so the editor draws it with powder conventions while
    we correctly call it unknown.

    While this app emitted an empty `##DATA TYPE=`, the editor's own
    `if (dataType)` guard sent it to PLAIN and the collision never fired.
    Preserving the string is right -- see the module docstring -- but it makes
    that collision live, and a single-crystal dataset is a reflection list, not
    a 1D diffractogram, so a powder rendering of it is meaningless rather than
    merely imprecise.

    The editor needs an entry ordered *ahead* of its `X-RAY DIFFRACTION` check
    before this ships. All 34 values of the converter app's dropdown
    (`converter_app/options.py`) were checked; this is the only collision.

    The assertion is deliberately about what we emit, since the editor's
    behaviour cannot be asserted from this repo.
    """
    core = JcampBaseConverter(
        write_with_datatype(tmp_path, 'SINGLE CRYSTAL X-RAY DIFFRACTION')
    )
    assert core.typ == '', 'we must not claim this is a known technique'
    assert core.datatype == 'SINGLE CRYSTAL X-RAY DIFFRACTION'
    assert 'X-RAY DIFFRACTION' in core.datatype, (
        'the substring the editor keys on; if this ever stops being true the '
        'cross-repo hazard is gone and this test can go with it'
    )
