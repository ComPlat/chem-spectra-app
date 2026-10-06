"""The block-aware reader, against real files.

Skipped on the pin this app currently ships, which has no `read_blocks`. Run
them against a venv built on ComPlat/nmrglue@b0802a0:

    <venv>/bin/python -m pytest tests/lib/converter/jcamp/test_reader.py

The point of these is that the reader reports what the *file* says, so they
assert block counts and order against files whose structure is known by
inspection rather than against another reader's output.
"""

import nmrglue as ng
import pytest

pytestmark = pytest.mark.skipif(
    not hasattr(ng.jcampdx, 'read_blocks'),
    reason='needs ComPlat/nmrglue@b0802a0 or later; the shipped pin has no read_blocks',
)

from chem_spectra.lib.converter.jcamp.reader import read_jcamp  # noqa: E402

NMR_LINK = './tests/fixtures/source/1H.dx'
SINGLE_BLOCK = './tests/fixtures/source/IR.dx'
MS = './tests/fixtures/source/ms/MS_ESI.jdx'


def test_an_nmr_link_file_reports_its_three_blocks_in_file_order():
    """`1H.dx` is a LINK wrapper, an FID block, then the spectrum. The old
    flat read collapsed these into one dict; the order here is the file's."""
    jcamp = read_jcamp(NMR_LINK)
    assert len(jcamp) == 3
    assert jcamp.datatypes == ['LINK', 'NMR FID', 'NMR SPECTRUM']


def test_a_single_block_file_reports_one_block():
    jcamp = read_jcamp(SINGLE_BLOCK)
    assert len(jcamp) == 1
    assert jcamp.datatypes == ['INFRARED SPECTRUM']


def test_the_measurement_block_is_chosen_by_file_order():
    """#291's rule -- first recognised datatype wins -- applied to blocks.

    On this file that is the third block, where the old `target_idx` was 1:
    the index into a merged list and the index of a block are different
    numbers, which is the whole reason for this layer.
    """
    jcamp = read_jcamp(NMR_LINK)
    target = jcamp.first_matching(lambda b: b.datatype == 'NMR SPECTRUM')
    assert target.index == 2


def test_ldr_reads_from_the_block_not_from_the_file():
    """Both NMR blocks carry `##UNITS=`; each must report its own.

    This is the defect the whole migration exists to fix: flattened, these two
    became one list and every caller had to guess an index into it.
    """
    jcamp = read_jcamp(NMR_LINK)
    fid, spectrum = jcamp[1], jcamp[2]
    assert fid.ldr('UNITS').split(',')[0].strip() == 'SECONDS'
    assert spectrum.ldr('UNITS').split(',')[0].strip() == 'HZ'


def test_ldrs_returns_every_occurrence_within_one_block():
    jcamp = read_jcamp(NMR_LINK)
    assert jcamp[2].ldrs('UNITS') == [jcamp[2].ldr('UNITS')]
    assert jcamp[0].ldrs('NOSUCHKEY') == []


def test_a_missing_ldr_is_none_rather_than_an_error():
    assert read_jcamp(SINGLE_BLOCK)[0].ldr('NOSUCHKEY') is None


@pytest.mark.parametrize('path, index, kind', [
    (SINGLE_BLOCK, 0, 'array'),     # (X++(Y..Y))
    (NMR_LINK, 2, 'ntuples'),       # NTUPLES
    (MS, 0, 'pairs'),               # (XY..XY) / PEAK TABLE
])
def test_data_comes_back_in_the_shape_the_block_declares(path, index, kind):
    block = read_jcamp(path)[index]
    data = block.data
    if kind == 'ntuples':
        assert isinstance(data, dict) and data.get('real')
    elif kind == 'pairs':
        assert data.shape[0] == 1 and data.shape[2] == 2
    else:
        assert data.ndim == 1 and len(data) > 0


def test_data_is_read_once_and_cached():
    """`getdataarray` re-parses the table on every call; callers touch `.data`
    in several places per file."""
    block = read_jcamp(SINGLE_BLOCK)[0]
    assert block.data is block.data


def test_a_block_with_no_data_reports_none():
    """The LINK wrapper carries LDRs but no table of its own."""
    assert read_jcamp(NMR_LINK)[0].data is None


def test_flat_ldrs_concatenates_every_block_in_file_order():
    """The metadata dump writes the whole file back as `###KEY= v1, v2, ...`,
    and nineteen golden files compare it byte for byte -- so the order and the
    multiplicity are both part of the contract.

    `1H.dx` has three blocks each declaring `##TITLE=`, which the dump emits as
    one joined value.
    """
    jcamp = read_jcamp(NMR_LINK)
    flat = jcamp.flat_ldrs()
    assert len(flat['TITLE']) == 3
    assert flat['TITLE'] == [b.ldr('TITLE') for b in jcamp]
    # both blocks' UNITS survive, in file order
    assert flat['UNITS'] == [jcamp[1].ldr('UNITS'), jcamp[2].ldr('UNITS')]


def test_flat_ldrs_drops_the_readers_own_bookkeeping():
    """`_parent` and friends are nmrglue's, not the file's, and must not reach
    the dump."""
    flat = read_jcamp(NMR_LINK).flat_ldrs()
    assert not [k for k in flat if k.startswith('_')]


def test_datatype_is_upper_cased_like_the_classifier_compares_it():
    assert read_jcamp(SINGLE_BLOCK)[0].datatype == 'INFRARED SPECTRUM'


# - - - a coordinate table may carry more than two columns - - -

def test_a_three_column_table_gives_its_x_and_y():
    """`(XYW..XYW)` keeps its third column since nmrglue `0aa0aa7`; before
    that the width was dropped and every table arrived two wide.

    `__block_pairs` asked for exactly two columns, so a three-column block
    fell through to `return data` and would have handed a three-dimensional
    array to the y series. It asks for two *or more* now and reads the first
    two, which is the same answer for an `(XY..XY)` table and the right one
    for `(XYW..XYW)`.
    """
    import numpy as np

    from chem_spectra.lib.converter.jcamp.technique import (
        JcampTechniqueConverter,
    )

    pairs = JcampTechniqueConverter._JcampTechniqueConverter__block_pairs
    xy = np.array([[[1.0, 10.0], [2.0, 20.0]]])          # (1, 2, 2)
    xyw = np.array([[[1.0, 10.0, 0.5], [2.0, 20.0, 0.6]]])  # (1, 2, 3)
    for data in (xy, xyw):
        x, y = pairs(data)
        assert list(x) == [1.0, 2.0]
        assert list(y) == [10.0, 20.0]
    assert pairs(np.array([[1.0, 10.0, 0.5], [2.0, 20.0, 0.6]])) is not None


def test_the_peak_table_width_does_not_disturb_x_and_y():
    """The same claim against the file that actually carries one.

    `CHI-224_10.jdx`'s NMR peak table reads `(1, 55, 3)` on this pin and read
    `(1, 55, 2)` before it. It is not that file's target block -- no fixture
    here has a three-column target -- so this pins the columns rather than any
    composed output.
    """
    from chem_spectra.lib.converter.jcamp.reader import read_jcamp
    from chem_spectra.lib.converter.jcamp.technique import (
        JcampTechniqueConverter,
    )

    blocks = read_jcamp('./tests/fixtures/source/CHI-224_10.jdx')
    table = blocks[2]
    assert table.datatype == 'NMR PEAK TABLE'
    assert table.data.shape[:2] == (1, 55)

    pairs = JcampTechniqueConverter._JcampTechniqueConverter__block_pairs
    x, y = pairs(table.data)
    assert len(x) == 55 and len(y) == 55
    assert x[0] == table.data[0][0][0]
    assert y[0] == table.data[0][0][1]
