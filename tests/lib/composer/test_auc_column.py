"""The AUC column in the integration table, and which techniques get it.

`integration_uses_auc_column` returns True either when the technique is HPLC
UV/VIS or when any integration item already carries `absoluteArea`. Every
integration fixture in the repo carries it, so the technique arm was dead
weight in the suite: flipping `auc_column` off in the registry failed
nothing at all.

These drive it with items that have no `absoluteArea`, which is the only
case where the technique decides. Per the plan's rule a zero sensitivity
count is a gap, not a pass.

Note `auc_column` is HPLC UV/VIS only, while `visual_split` also covers
plain UV/VIS -- a wider set. That asymmetry is the reason they are two
fields, and is asserted here.
"""

import re

from chem_spectra.lib.composer.technique import TechniqueComposer
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter

SOURCE = './tests/fixtures/source/hplc/chromatogram.jdx'

# no 'absoluteArea' key: the technique is what decides the column count
ITEMS_WITHOUT_ABSOLUTE_AREA = [
    {'xL': 1.0, 'xU': 2.0, 'area': 10.0},
    {'xL': 3.0, 'xU': 4.0, 'area': 20.0},
]


def _composer(datatype, tmp_path):
    body = open(SOURCE).read()
    body = body.replace(
        re.search(r'##DATA TYPE=.*', body).group(0), '##DATA TYPE=' + datatype, 1)
    target = tmp_path / 'probe.jdx'
    target.write_text(body)
    return TechniqueComposer(JcampTechniqueConverter(JcampBaseConverter(str(target))))


def _columns(line):
    return len(line.strip().strip('()').split(','))


def test_hplc_uv_vis_gets_the_auc_column(tmp_path):
    composer = _composer('HPLC UV/VIS SPECTRUM', tmp_path)
    assert composer._technique().auc_column is True
    lines = composer._build_integration_lines(ITEMS_WITHOUT_ABSOLUTE_AREA)
    assert lines and all(_columns(line) == 4 for line in lines)


def test_plain_uv_vis_does_not(tmp_path):
    """The narrower set: UV/VIS supports visual splits but not the AUC column."""
    composer = _composer('UV/VIS SPECTRUM', tmp_path)
    assert composer._technique().auc_column is False
    assert composer._technique().visual_split is True
    lines = composer._build_integration_lines(ITEMS_WITHOUT_ABSOLUTE_AREA)
    assert lines and all(_columns(line) == 3 for line in lines)


def test_an_item_carrying_absolute_area_gets_the_column_regardless(tmp_path):
    """The other arm of the or: the data can ask for it too."""
    composer = _composer('UV/VIS SPECTRUM', tmp_path)
    items = [dict(ITEMS_WITHOUT_ABSOLUTE_AREA[0], absoluteArea=5.0)]
    lines = composer._build_integration_lines(items)
    assert lines and _columns(lines[0]) == 4


def test_the_two_uv_vis_fields_cover_different_sets():
    from chem_spectra.lib.converter.jcamp.techniques import SPECTRUM_TECHNIQUES
    auc = {k for k, t in SPECTRUM_TECHNIQUES.items() if t.auc_column}
    split = {k for k, t in SPECTRUM_TECHNIQUES.items() if t.visual_split}
    assert auc == {'HPLC UVVIS'}
    assert split == {'HPLC UVVIS', 'UVVIS'}
    assert auc < split
