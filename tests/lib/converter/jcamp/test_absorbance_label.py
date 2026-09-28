"""How a y-axis declaring absorbance is labelled, per technique.

`__set_label` rewrites an absorbance y-axis to TRANSMITTANCE for every
technique except UV/VIS. These tests pin that, including the case that looks
wrong, so the behaviour cannot drift unnoticed while the question is open.

**Suspected defect, deliberately preserved.** 'HPLC UVVIS' is a separate
registry key from 'UVVIS', so an HPLC UV/VIS file declaring
##YUNITS=ABSORBANCE is relabelled TRANSMITTANCE -- the inverse quantity.
The shipped HPLC fixture declares mAU and so escapes it, but the path is
reachable. Whether HPLC UV/VIS should keep its absorbance label is a domain
call for whoever runs the instrument, not the refactor's to make, so it is
pinned rather than fixed. See CHANGELOG.refactor-finish-flag-migration.md.
"""

import re

import pytest

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter

SOURCE = './tests/fixtures/source/JPK-948.jdx'   # a real UV/VIS, YUNITS=ABSORBANCE


def _probe(datatype, tmp_path, yunits='ABSORBANCE'):
    body = open(SOURCE).read()
    body = body.replace(
        re.search(r'##DATA TYPE=.*', body).group(0), '##DATA TYPE=' + datatype, 1)
    body = body.replace(
        re.search(r'##YUNITS=.*', body).group(0), '##YUNITS=' + yunits, 1)
    target = tmp_path / 'probe.jdx'
    target.write_text(body)
    return JcampTechniqueConverter(JcampBaseConverter(str(target)))


def test_uv_vis_keeps_its_absorbance_label(tmp_path):
    converter = _probe('UV/VIS SPECTRUM', tmp_path)
    assert converter.technique.absorbance_label is True
    assert converter.label['y'] == 'ABSORBANCE'


@pytest.mark.parametrize('datatype', [
    'INFRARED SPECTRUM',
    'RAMAN SPECTRUM',
])
def test_other_techniques_are_rewritten_to_transmittance(datatype, tmp_path):
    converter = _probe(datatype, tmp_path)
    assert converter.technique.absorbance_label is False
    assert converter.label['y'] == 'TRANSMITTANCE'


def test_hplc_uv_vis_is_rewritten_too_which_is_the_open_question(tmp_path):
    """Pins a suspected defect rather than fixing it -- see the module docstring.

    If this is ever judged wrong, the fix is one field on the HPLC UVVIS
    registry entry, and this test is the one to invert.
    """
    converter = _probe('HPLC UV/VIS SPECTRUM', tmp_path)
    assert converter.technique.absorbance_label is False
    assert converter.label['y'] == 'TRANSMITTANCE'


def test_a_non_absorbance_label_is_left_alone(tmp_path):
    """The rewrite is keyed on the declared units, not on the technique."""
    converter = _probe('HPLC UV/VIS SPECTRUM', tmp_path, yunits='mAU')
    assert converter.label['y'] == 'mAU'
