"""Size-exclusion and differential-scanning outputs, neither of which had a test.

No fixture in the repo declares either technique, so the `is_sec` / `is_dsc`
branches in TechniqueComposer were never executed by the suite -- before this
refactor as well as after. All four descriptor fields they migrated to
reported a zero sensitivity count, which the plan treats as a gap rather than
a pass.

These drive both techniques with a relabelled probe carrying the LDRs each
one reads, and assert the four distinct outputs:

- `sec_headers`  -- ##MN/##MW/##MP/##D in the JCAMP metadata
- `dsc_metadata` -- melting point and Tg in the JCAMP metadata
- `info_box`     -- which annotation box is drawn on the plot, if any
"""

import re

from chem_spectra.lib.composer.technique import TechniqueComposer
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter

SOURCE = './tests/fixtures/source/hplc/chromatogram.jdx'

SEC_LDRS = '##MN=4200\n##MW=8600\n##MP=7100\n##D=2.05\n'
DSC_LDRS = '##MELTINGPOINT=131.4\n##TG=58.2\n'


def _composer(datatype, extra_ldrs, tmp_path, params=None):
    body = open(SOURCE).read()
    header = re.search(r'##DATA TYPE=.*', body).group(0)
    body = body.replace(header, '##DATA TYPE=' + datatype + '\n' + extra_ldrs.rstrip(), 1)
    target = tmp_path / 'probe.jdx'
    target.write_text(body)
    base = JcampBaseConverter(str(target), params) if params else JcampBaseConverter(str(target))
    return TechniqueComposer(JcampTechniqueConverter(base))


def _info_box_text(composer):
    """Capture whatever __generate_info_box writes, without rendering a file."""
    captured = {}

    class FakeAxes:
        transAxes = object()

        def text(self, *args, **kwargs):
            captured['text'] = args[2] if len(args) > 2 else kwargs.get('s')

    class FakePlot:
        def gca(self):
            return FakeAxes()

    composer._TechniqueComposer__generate_info_box(FakePlot())
    return captured.get('text')


# - - - sec_headers - - -

def test_sec_writes_its_header_block(tmp_path):
    composer = _composer('SIZE EXCLUSION CHROMATOGRAPHY', SEC_LDRS, tmp_path)
    assert composer._technique().sec_headers is True
    meta = ''.join(composer.meta)
    for ldr in ('##MN=4200', '##MW=8600', '##MP=7100', '##D=2.05'):
        assert ldr in meta


def test_a_non_sec_technique_writes_no_sec_block(tmp_path):
    composer = _composer('HPLC UV/VIS SPECTRUM', SEC_LDRS, tmp_path)
    assert composer._technique().sec_headers is False
    meta = ''.join(composer.meta)
    assert '##MN=4200' not in meta


# - - - info_box - - -

def test_sec_draws_the_molecular_weight_box(tmp_path):
    composer = _composer('SIZE EXCLUSION CHROMATOGRAPHY', SEC_LDRS, tmp_path)
    assert composer._technique().info_box == 'sec'
    text = _info_box_text(composer)
    assert text is not None
    assert 'MN=4200' in text and 'D=2.05' in text
    assert 'MELTING POINT' not in text


def test_dsc_draws_the_melting_point_box(tmp_path):
    composer = _composer('DIFFERENTIAL SCANNING CALORIMETRY', DSC_LDRS, tmp_path)
    assert composer._technique().info_box == 'dsc'
    text = _info_box_text(composer)
    assert text is not None
    assert 'MELTING POINT=131.4' in text and 'TG=58.2' in text
    assert 'MN=' not in text


def test_every_other_technique_draws_no_box(tmp_path):
    composer = _composer('HPLC UV/VIS SPECTRUM', SEC_LDRS, tmp_path)
    assert composer._technique().info_box == ''
    assert _info_box_text(composer) is None


# - - - dsc_metadata - - -

def test_dsc_writes_melting_point_and_tg_from_the_file(tmp_path):
    composer = _composer('DIFFERENTIAL SCANNING CALORIMETRY', DSC_LDRS, tmp_path)
    assert composer._technique().dsc_metadata is True
    meta = ''.join(composer.meta)
    assert '131.4' in meta and '58.2' in meta


def test_a_non_dsc_technique_writes_neither(tmp_path):
    composer = _composer('HPLC UV/VIS SPECTRUM', DSC_LDRS, tmp_path)
    assert composer._technique().dsc_metadata is False
