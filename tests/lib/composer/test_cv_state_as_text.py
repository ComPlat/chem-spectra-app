"""A CV state handed to the composer as JSON text is read, not lost.

`parse_params` hands the converter a dict, but the composer accepts the
state as text too -- and `composer/technique.py` called `json.loads` without
importing `json`. In `__cv_display_info` the NameError was swallowed and the
y scale silently dropped; in the CV metadata it raised.
"""

import json

from chem_spectra.lib.composer.technique import TechniqueComposer
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter

SOURCE = './tests/fixtures/source/cyclicvoltammetry/RCV_LSH-R444_full+Fc.jdx'

STATE = {
    'spectraList': [{'list': [], 'shift': {}, 'hasRefPeak': False}],
    'cvDisplay': {'yScaleFactor': 2},
    'areaValue': '0.07',
    'areaUnit': 'cm²',
    'useCurrentDensity': True,
}


def _composer_with_text_state():
    converter = JcampTechniqueConverter(JcampBaseConverter(SOURCE))
    converter.params['cyclicvolta'] = json.dumps(STATE)
    return TechniqueComposer(converter)


def test_the_cv_metadata_reads_a_text_state():
    meta = ''.join(_composer_with_text_state().meta)
    assert '##$CSWEAREAVALUE=0.07\n' in meta
    assert '##$CSCURRENTMODE=DENSITY\n' in meta


def test_the_display_scale_is_read_from_a_text_state():
    composer = _composer_with_text_state()
    info = composer._TechniqueComposer__cv_display_info()
    assert info['cv_display'] == {'yScaleFactor': 2}
