"""Every mass-spectrum run in a file, whatever each block calls itself.

An MS file can hold one run per block. The flat read merged every block's
pages; reading per block, `make_ms_data_xsys` gathers them on purpose. It
gathered them by the target's exact `##DATA TYPE=` string, so a file mixing
`MASS SPECTRUM` and `CONTINUOUS MASS SPECTRUM` -- both mapped to MS -- lost the
other spelling's runs, silently.
"""

import re

import nmrglue as ng
import pytest

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ms import JcampMSConverter

pytestmark = pytest.mark.skipif(
    not hasattr(ng.jcampdx, 'read_blocks'),
    reason='needs ComPlat/nmrglue@b0802a0 or later; the shipped pin has no read_blocks',
)

SOURCE = './tests/fixtures/source/ms/svs813f1.jdx'   # one block, 24 runs


def _two_blocks(tmp_path, second_datatype):
    text = open(SOURCE).read()
    block = text[text.index('##TITLE='):]
    block = block[:block.index('##END=')] + '##END=\n'
    second = re.sub(r'##DATA TYPE=.*', '##DATA TYPE=' + second_datatype,
                    block, count=1)
    target = tmp_path / 'two_runs.jdx'
    target.write_text('##TITLE=two runs\n##JCAMP-DX=5.00\n##DATA TYPE=LINK\n'
                      '##BLOCKS=2\n' + block + second + '##END=\n')
    return target


@pytest.mark.parametrize('second', ['MASS SPECTRUM',
                                    'CONTINUOUS MASS SPECTRUM'])
def test_every_ms_block_contributes_its_runs(tmp_path, second):
    path = _two_blocks(tmp_path, second)
    converter = JcampMSConverter(JcampBaseConverter(str(path), {}))
    assert len(converter.data) == 48


def test_a_tic_is_not_gathered_as_a_mass_spectrum(tmp_path):
    """It maps to LC/MS, not MS, so the target's runs stay its own."""
    path = _two_blocks(tmp_path, 'TOTAL ION CHROMATOGRAM')
    converter = JcampMSConverter(JcampBaseConverter(str(path), {}))
    assert len(converter.data) == 24
