"""Avantes `AvaSoft` exports are UV/VIS, not an unrecognised datatype.

Avantes spectrometer software writes `##DATA TYPE= AvaSoft` -- the product
name rather than a measurement name -- so classification did not recognise it.

Before #291 that raised. #291 stopped the crash by routing unrecognised
datatypes to the generic curve, which is why the file renders today, but it
still arrived unclassified: `typ=''`, no technique, and the generic peak
threshold of 0.5 instead of UV/VIS's 0.05.

PR #243 proposed this alias in July 2025 and was closed unmerged. It mapped
AvaSoft to **HPLC UVVIS**, which is wrong: the file carries
`##XUNITS= NANOMETERS`, a wavelength scan, where an HPLC chromatogram is a
time series (`##XUNITS=MINUTES`). Under the HPLC key it would also inherit
that technique's AUC-column and visual-split integration behaviour.

Reported as issue #242.
"""

import json
import os

import pytest

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter

DATA_TYPE_JSON = os.path.join(
    os.path.dirname(
        __import__('chem_spectra.lib.converter.jcamp.base',
                   fromlist=['base']).__file__),
    'data_type.json',
)


def _mapping():
    with open(DATA_TYPE_JSON) as handle:
        return json.load(handle)['datatypes']


def test_avasoft_is_a_uv_vis_alias():
    assert 'AvaSoft' in _mapping()['UVVIS']


def test_avasoft_is_not_mapped_to_hplc():
    """The distinction PR #243 got wrong: nanometres, not minutes."""
    assert 'AvaSoft' not in _mapping()['HPLC UVVIS']


@pytest.fixture
def avasoft_file(tmp_path):
    """Minimal reproduction of the Avantes header from issue #242."""
    xs = [400.0 + i for i in range(50)]
    ys = [0.1] * 50
    ys[20] = 0.9
    body = [
        '##TITLE= Vial 35\n',
        '##JCAMP-DX= 4.24\n',
        '##DATA TYPE= AvaSoft\n',
        '##ORIGIN=AvaSpec 7213289SP\n',
        '##OWNER=\n',
        '##NPOINTS= {}\n'.format(len(xs)),
        '##XUNITS= NANOMETERS\n',
        '##YUNITS= TRANSMITTANCE\n',
        '##FIRSTX= {}\n'.format(xs[0]),
        '##LASTX= {}\n'.format(xs[-1]),
        '##XFACTOR= 1.0\n',
        '##YFACTOR= 1.0\n',
        '##XYPOINTS= (XY..XY)\n',
    ]
    body += [' {:.2f},{:.5f}\n'.format(x, y) for x, y in zip(xs, ys)]
    body.append('##END=\n')
    target = tmp_path / 'avasoft.dx'
    target.write_text(''.join(body))
    return str(target)


def test_an_avasoft_file_classifies_as_uv_vis(avasoft_file):
    base = JcampBaseConverter(avasoft_file)
    assert base.typ == 'UVVIS'
    assert base.technique.key == 'UVVIS'


def test_it_gets_the_uv_vis_threshold_not_the_generic_one(avasoft_file):
    """The observable consequence of being unclassified: 0.5 rather than 0.05."""
    converter = JcampTechniqueConverter(JcampBaseConverter(avasoft_file))
    assert converter.threshold == 0.05
