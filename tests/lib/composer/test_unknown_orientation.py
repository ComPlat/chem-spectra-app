"""An unrecognised file is drawn ascending, and the registry must say so.

`UNKNOWN_TECHNIQUE` used to inherit `x_reversed=True` from the
`SpectrumTechnique` default. That was never a decision about unrecognised
files -- every one of the eighteen mapped entries passes the flag
explicitly, so the default had exactly one reader, and what it said was
"draw an unidentified x quantity the way NMR draws ppm".

The cost is cross-repo. react-spectra-editor#336 draws its PLAIN layout
ascending (`Format.isNonReversedXLayout`), and PLAIN is where every
unrecognised `##DATA TYPE=` lands. Left reversed here, the preview image
this app renders and the curve the editor draws for the same file are
mirror images of each other.

Measured rather than restated, following `test_ms_orientation.py`: the
registry value is only worth what the render does with it. The axis limits
have to be read *during* `savefig`, since `tf_img` calls `plt.clf()`
immediately afterwards.

The probe file is `IR.dx` with only its `##DATA TYPE=` swapped for `SQUID`
-- a real datatype chemotion-converter-app emits and this app deliberately
does not map. Its points run 3997 -> 374, i.e. descending in storage, and
`SQUID` is not an `em_wave` technique so nothing reorders them. So the
assertion below is a statement about the drawn axis alone, which is the
thing the editor has to agree with.
"""

import chem_spectra.lib.composer.technique as technique_module
from chem_spectra.lib.composer.technique import TechniqueComposer
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter
from chem_spectra.lib.converter.jcamp.techniques import (
    SPECTRUM_TECHNIQUES, UNKNOWN_TECHNIQUE,
)
from tests.lib.converter.jcamp.test_unknown_datatype_preserved import (
    write_with_datatype,
)


def _rendered_xlim(composer):
    captured = {}
    real_savefig = technique_module.plt.savefig

    def spy_savefig(*args, **kwargs):
        captured['xlim'] = technique_module.plt.gca().get_xlim()
        return real_savefig(*args, **kwargs)

    technique_module.plt.savefig = spy_savefig
    try:
        composer.tf_img().close()
    finally:
        technique_module.plt.savefig = real_savefig
    return captured['xlim']


def _compose(path):
    return TechniqueComposer(JcampTechniqueConverter(JcampBaseConverter(path)))


def test_the_registry_entry_is_not_reversed():
    assert UNKNOWN_TECHNIQUE.x_reversed is False


def test_an_unrecognised_file_renders_ascending(tmp_path):
    """The measurement the registry value has to agree with."""
    path = write_with_datatype(tmp_path, 'SQUID')
    low, high = _rendered_xlim(_compose(path))
    assert low < high, 'an unrecognised file rendered descending'


def test_the_stored_points_still_descend(tmp_path):
    """Why the test above proves something.

    The probe's own data runs high -> low. If a future change reordered it,
    the assertion would start passing for the wrong reason.
    """
    path = write_with_datatype(tmp_path, 'SQUID')
    core = JcampTechniqueConverter(JcampBaseConverter(path))
    assert core.technique is UNKNOWN_TECHNIQUE
    assert core.xs[0] > core.xs[-1]


def test_a_recognised_reversed_technique_is_unaffected():
    """The flip must not reach anything that chose its orientation.

    IR is the neighbouring case: same fixture, same generic curve styling,
    and it stays reversed because its entry says so.
    """
    assert SPECTRUM_TECHNIQUES['INFRARED'].x_reversed is True
    low, high = _rendered_xlim(_compose('./tests/fixtures/source/IR.dx'))
    assert low > high, 'an infrared spectrum rendered ascending'
