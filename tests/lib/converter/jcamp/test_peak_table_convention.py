"""The two peak tables are told apart by position, and nothing checks it.

A ChemSpectra output carries its peaks in two `##DATA CLASS=PEAKTABLE` blocks.
Each one declares which it is -- `##$CSCATEGORY=EDIT_PEAK` or `AUTO_PEAK` --
but the reader ignores that and goes by position: `PEAKTABLE[0]` is the edit
table and `PEAKTABLE[1]` is the auto one, in the order the composer wrote them
(`jcamp/technique.py`, `__read_edit_peaks` / `__read_auto_peaks`).

That works only because nmrglue's current flat read happens to accumulate the
blocks in file order. The nmrglue migration replaces that read, so the
assumption stops being free. These tests pin the convention to the categories
the file actually declares, so a migration that reorders or drops a block
fails here rather than silently swapping a chemist's hand-picked peaks for the
automatic ones.

`1H.edit.jdx` is the only fixture with both tables and no test used it.
"""

import re

import pytest

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter

SOURCE = './tests/fixtures/source/1H.edit.jdx'


def peaks_declared_as(category):
    """Read the peak table that names itself `category`, straight from the file.

    Deliberately independent of the app's reader: the point is to compare the
    positional convention against what the file declares, so this parse must
    not share the assumption under test.
    """
    with open(SOURCE) as handle:
        text = handle.read()

    blocks = text.split('##DATA TYPE=NMRPEAKTABLE')
    for block in blocks[1:]:
        if f'##$CSCATEGORY={category}' not in block:
            continue
        body = block.split('##PEAKTABLE= (XY..XY)')[1]
        pairs = []
        for line in body.splitlines():
            line = line.strip()
            if not line or line.startswith('##') or line.startswith('$$'):
                if line.startswith('##END'):
                    break
                continue
            if not re.match(r'^-?[\d.]+\s*,', line):
                continue
            x, y = line.split(',')[:2]
            pairs.append((float(x), float(y)))
        return pairs
    raise AssertionError(f'no {category} block in {SOURCE}')


@pytest.fixture
def converter():
    return JcampTechniqueConverter(JcampBaseConverter(SOURCE))


def test_the_fixture_really_declares_both_tables():
    """Guards the two tests below: if the fixture changes, they must not
    silently start comparing empty lists."""
    assert len(peaks_declared_as('EDIT_PEAK')) > 0
    assert len(peaks_declared_as('AUTO_PEAK')) > 0
    assert peaks_declared_as('EDIT_PEAK') != peaks_declared_as('AUTO_PEAK')


def test_the_first_peak_table_is_the_edit_table(converter):
    declared = peaks_declared_as('EDIT_PEAK')
    assert converter.edit_peaks is not None, 'edit peaks were not read at all'
    assert converter.edit_peaks['x'] == pytest.approx([x for x, _ in declared])
    assert converter.edit_peaks['y'] == pytest.approx([y for _, y in declared])


def test_the_second_peak_table_is_the_auto_table(converter):
    declared = peaks_declared_as('AUTO_PEAK')
    assert converter.auto_peaks is not None, 'auto peaks were not read at all'
    assert converter.auto_peaks['x'] == pytest.approx([x for x, _ in declared])
    assert converter.auto_peaks['y'] == pytest.approx([y for _, y in declared])


def test_the_two_tables_are_not_the_same_one(converter):
    """The failure mode worth catching: both handles pointing at one block."""
    assert converter.edit_peaks['x'] != converter.auto_peaks['x']
