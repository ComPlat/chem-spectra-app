"""The header reads in JcampTechniqueConverter fall back on what they expect, and
on nothing else.

Each of these reads used a bare `except:`, which turned any failure into a
default -- a missing record, but equally a bug. They now catch the errors a
missing or malformed record raises, so those still fall back, and anything
else surfaces.

The files are the repo's own fixtures with one header edited, written to
tmp_path.
"""
import pytest

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter

SOURCE = './tests/fixtures/source/'


def _edited(tmp_path, name, old, new, count=-1):
    with open(SOURCE + name, encoding='utf-8', errors='ignore') as f:
        text = f.read()
    assert old in text
    path = tmp_path / name.split('/')[-1]
    path.write_text(text.replace(old, new, count), encoding='utf-8')
    return str(path)


def _convert(path):
    return JcampTechniqueConverter(JcampBaseConverter(path))


def test_an_unreadable_block_count_counts_as_one(tmp_path):
    converter = _convert(_edited(tmp_path, '1H.dx', '##BLOCKS=   2', '##BLOCKS= two'))

    assert converter.block_count == 1


def test_an_unreadable_observe_frequency_falls_back_to_sfo1(tmp_path):
    path = _edited(tmp_path, '1H.dx',
                   '##.OBSERVE FREQUENCY= 400.13240078',
                   '##.OBSERVE FREQUENCY= unknown')

    converter = _convert(path)

    assert converter.obs_freq == pytest.approx(400.13240078)
    assert converter.ys is not None


def test_an_unreadable_xfactor_falls_back_to_the_default(tmp_path):
    converter = _convert(_edited(tmp_path, 'IR.dx', '##XFACTOR=1', '##XFACTOR=one'))

    assert converter.factor['x'] == 1.0


def test_a_malformed_peak_row_drops_that_table_and_still_converts(tmp_path):
    """A stored table with a row that is not a number is dropped, as before,
    rather than raising; the peaks are then picked from the trace."""
    stored = [5.5, 6.6]  # the file's own edit table
    whole = _convert(SOURCE + 'edit/edit_1H.dx')
    assert whole.edit_peaks['x'] == stored

    path = _edited(tmp_path, 'edit/edit_1H.dx',
                   '##PEAK ASSIGNMENTS=(XYA)\n',
                   '##PEAK ASSIGNMENTS=(XYA)\n(not, a, row)\n')
    cut = _convert(path)

    assert cut.edit_peaks['x'] != stored
    assert cut.auto_peaks
    assert cut.ys is not None


class _BrokenBlocks(dict):
    def __getitem__(self, key):
        if key == 'BLOCKS':
            raise RuntimeError('a bug, not a missing record')
        return super().__getitem__(key)


def test_an_unexpected_error_in_a_header_read_is_not_swallowed():
    """The block count falls back when BLOCKS is absent (KeyError) or
    unreadable (ValueError). A different error there is a bug, and it now
    surfaces instead of becoming a block count of one."""
    base = JcampBaseConverter(SOURCE + 'IR.dx')
    base.dic = _BrokenBlocks(base.dic)

    with pytest.raises(RuntimeError, match='a bug'):
        JcampTechniqueConverter(base)
