"""The stored point order follows the descriptor, not the header shape.

`JcampTechniqueConverter.__read_xs` reverses both axes when a file stores x
the other way round from how the technique is drawn. Two questions decide
it, and they used to be answered by one field:

* **which way is it drawn** -- `x_reversed`, the same field the preview
  image uses and, since react-spectra-editor#336, the editor's own axis list;
* **may the stored order be changed at all** -- `store_in_drawn_order`,
  opt-in, because for some techniques the traversal order is data. A
  sorption isotherm's branch is labelled ADSORPTION or DESORPTION by which
  way x runs (`tests/model/test_sorption_branches.py`), and a cyclic
  voltammogram's sweep likewise.

Both used to come from `em_wave`, which is a header-shape grouping -- its
other reader decides whether to emit the NMR header block -- and says
nothing about an axis. The provenance shows it was never meant to:
`40da642` (2019) wrote the reversal as `self.typ == 'INFRARED'`; `af68fdd`
widened it to the grouping when Raman arrived; `387083a` added UV-Vis to
the grouping while adding the `UV/VIS SPECTRUM` identifier, a threshold and
a bin count, with nothing about direction. UV-Vis ascends by convention, so
it was stored backwards.

Scope: this code is on the `(X++(Y..Y))` path only. `__read_xs` returns from
`make_ni_data_xs` for `(XY..XY)` files, which is every cyclic-voltammetry
and HPLC fixture -- and `JPK-948.jdx`, the repository's only real UV-Vis
file, which already ascends because the instrument wrote it that way.

The probe is `IR.dx` with its datatype swapped, and optionally its declared
ends exchanged so the same samples read the other way. It must not be an
`(XY..XY)` file: `__read_xs` returns early for those and never reaches this.
"""

import re

import pytest

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter

SOURCE = './tests/fixtures/source/IR.dx'   # stores x descending, 3997 -> 374


def _probe(datatype, tmp_path, ascending=False, name='probe.dx'):
    body = open(SOURCE).read()
    header = re.search(r'##DATA TYPE=.*', body).group(0)
    body = body.replace(header, '##DATA TYPE=' + datatype, 1)
    if ascending:
        body = body.replace('##FIRSTX=3997.453', '##FIRSTX=@@') \
                   .replace('##LASTX=373.96442', '##LASTX=3997.453') \
                   .replace('##FIRSTX=@@', '##FIRSTX=373.96442') \
                   .replace('##DELTAX=-1.4165319', '##DELTAX=1.4165319')
    target = tmp_path / name
    target.write_text(body)
    return JcampTechniqueConverter(JcampBaseConverter(str(target)))


# - - - techniques that opt in - - -

@pytest.mark.parametrize('datatype', [
    'INFRARED SPECTRUM',
    'RAMAN SPECTRUM',
])
def test_a_reversed_technique_is_stored_descending(datatype, tmp_path):
    """The case the reversal was written for, from an ascending file."""
    converter = _probe(datatype, tmp_path, ascending=True)
    assert converter.technique.store_in_drawn_order is True
    assert converter.xs[0] > converter.xs[-1], (
        '%s must be stored descending' % datatype)


@pytest.mark.parametrize('ascending', [False, True])
def test_uvvis_is_stored_ascending(ascending, tmp_path):
    """What this changes: UV-Vis ascends by convention.

    Under the `em_wave` rule the descending source was left descending --
    that rule only ever forced one direction -- and the ascending one was
    turned round into it.
    """
    converter = _probe('UV/VIS SPECTRUM', tmp_path, ascending=ascending)
    assert converter.technique.x_reversed is False
    assert converter.xs[0] < converter.xs[-1]


def test_both_axes_move_together(tmp_path):
    """Or the spectrum is silently mirrored.

    Exchanging the declared ends without reversing `ys` would leave every
    intensity at the wrong wavelength, and the assertions above would still
    pass.
    """
    turned = _probe('UV/VIS SPECTRUM', tmp_path, name='uvvis.dx')
    left = _probe('HPLC UV/VIS SPECTRUM', tmp_path, name='hplc.dx')
    assert list(turned.ys) == list(reversed(list(left.ys)))
    assert turned.xs[0] == pytest.approx(left.xs[-1])
    assert turned.xs[-1] == pytest.approx(left.xs[0])


# - - - techniques that do not - - -

def test_a_technique_outside_the_rule_is_left_alone(tmp_path):
    """The control: the same bytes, a technique that has not opted in."""
    converter = _probe('HPLC UV/VIS SPECTRUM', tmp_path, ascending=True)
    assert converter.technique.store_in_drawn_order is False
    assert converter.xs[0] < converter.xs[-1]


def test_a_sorption_trace_keeps_the_direction_it_was_given(tmp_path):
    """Why the rule is opt-in rather than `x_reversed` alone.

    `SORPTION-DESORPTION MEASUREMENT` is drawn ascending, so a descending
    trace would be turned round -- and `tf_combine` reads that direction to
    label the branch ADSORPTION or DESORPTION. Normalising it relabels every
    desorption branch as its opposite.
    """
    converter = _probe('SORPTION-DESORPTION MEASUREMENT', tmp_path)
    assert converter.technique.x_reversed is False
    assert converter.technique.store_in_drawn_order is False
    assert converter.xs[0] > converter.xs[-1]


@pytest.mark.parametrize('fixture', [
    './tests/fixtures/source/1H.dx',
    './tests/fixtures/source/13C-CPD.dx',
    './tests/fixtures/source/mnova/MNOVA_SVS_13C.jdx',
])
def test_nmr_keeps_its_order(fixture):
    """NMR has not opted in either, and its fixtures already descend."""
    converter = JcampTechniqueConverter(JcampBaseConverter(fixture))
    assert converter.xs[0] > converter.xs[-1]


def test_the_xy_pair_path_is_untouched():
    """`JPK-948.jdx` is `(XY..XY)`, so `__read_xs` returns before any of
    this. It ascends because the instrument wrote it that way -- the order
    an equidistant UV-Vis file now also gets."""
    base = JcampBaseConverter('./tests/fixtures/source/JPK-948.jdx')
    assert base.data_format == '(XY..XY)'
    converter = JcampTechniqueConverter(base)
    assert converter.xs[0] < converter.xs[-1]
