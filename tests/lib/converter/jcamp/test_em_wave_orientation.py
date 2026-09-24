"""An em-wave spectrum stored ascending is flipped to descending.

`JcampTechniqueConverter.__read_xs` reverses both axes when the technique is
an electromagnetic-wave one and the file stores x ascending -- infrared,
Raman and UV/VIS are conventionally plotted with wavenumber decreasing.

No fixture exercised this. Every IR fixture in the repo already stores
FIRSTX > LASTX, so the branch never fired, and flipping `em_wave` off in the
registry failed nothing at all. These tests drive it with an ascending file,
so the field is load-bearing rather than merely present.

The probe is `IR.dx` with FIRSTX and LASTX exchanged -- declaring the same
samples run the other way. Note it must not be an `(XY..XY)` file:
`__read_xs` returns early for those and never reaches the flip.
"""

import re

import pytest

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter

SOURCE = './tests/fixtures/source/IR.dx'


def _ascending_probe(datatype, tmp_path):
    body = open(SOURCE).read()
    header = re.search(r'##DATA TYPE=.*', body).group(0)
    body = body.replace(header, '##DATA TYPE=' + datatype, 1)
    # exchange the declared ends so the file reads ascending
    body = body.replace('##FIRSTX=3997.453', '##FIRSTX=@@') \
               .replace('##LASTX=373.96442', '##LASTX=3997.453') \
               .replace('##FIRSTX=@@', '##FIRSTX=373.96442') \
               .replace('##DELTAX=-1.4165319', '##DELTAX=1.4165319')
    target = tmp_path / 'probe.dx'
    target.write_text(body)
    return JcampTechniqueConverter(JcampBaseConverter(str(target)))


@pytest.mark.parametrize('datatype', [
    'INFRARED SPECTRUM',
    'RAMAN SPECTRUM',
    'UV/VIS SPECTRUM',
])
def test_em_wave_ascending_file_is_flipped(datatype, tmp_path):
    converter = _ascending_probe(datatype, tmp_path)
    assert converter.technique.em_wave is True
    assert converter.xs[0] > converter.xs[-1], (
        '%s stored ascending must render descending' % datatype)


def test_non_em_wave_ascending_file_is_left_alone(tmp_path):
    """The control: the same bytes, a technique that is not em-wave."""
    converter = _ascending_probe('HPLC UV/VIS SPECTRUM', tmp_path)
    assert converter.technique.em_wave is False
    assert converter.xs[0] < converter.xs[-1]


def test_em_wave_flip_reverses_the_y_series_too(tmp_path):
    """Both axes move together, or the spectrum is silently mirrored."""
    em_wave = _ascending_probe('INFRARED SPECTRUM', tmp_path)
    plain = _ascending_probe('HPLC UV/VIS SPECTRUM', tmp_path)
    assert list(em_wave.ys) == list(reversed(list(plain.ys)))
