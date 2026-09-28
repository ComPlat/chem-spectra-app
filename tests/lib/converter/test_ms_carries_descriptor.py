"""Every MS converter carries the MS descriptor.

Three converters feed `MSComposer`, and all three state `typ = 'MS'`:

- `JcampMSConverter`, from its base converter's classification;
- `CdfMSConverter`, from `CdfBaseConverter`, which hardcodes it;
- `MSConverter` (RAW / mzML / mzXML), which hardcodes it literally.

The last two never go through JCAMP classification at all, which is what
made "which technique is this?" look like an open question for them. It is
not: `typ` is already stated, so `technique_for` answers it.

Carrying the descriptor is what lets `BaseComposer._technique()` resolve
these cores directly. Before this, `JcampMSConverter` copied `non_nmr` and
the other two set nothing, so all three fell through to the
`getattr(core, 'non_nmr', True)` default and resolved to
`UNKNOWN_TECHNIQUE`.

**Coverage gap, stated rather than papered over:** `CdfMSConverter` is not
asserted here. No `.cdf` fixture exists and `test_to_composer_cdf` in
`tests/model/test_transformer.py` is an empty stub, so nothing in the suite
constructs it. Its descriptor assignment is unverified by execution.
"""

from werkzeug.datastructures import FileStorage

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ms import JcampMSConverter
from chem_spectra.lib.converter.jcamp.techniques import SPECTRUM_TECHNIQUES
from chem_spectra.lib.converter.ms import MSConverter
from chem_spectra.lib.composer.ms import MSComposer


def test_jcamp_ms_converter_carries_the_descriptor():
    core = JcampMSConverter(JcampBaseConverter('./tests/fixtures/source/ms/svs813f1.jdx'))
    assert core.technique is SPECTRUM_TECHNIQUES['MS']


def test_mzml_converter_carries_the_descriptor():
    """The path that never sees data_type.json."""
    with open('./tests/fixtures/source/ms/svs813f1.mzML', 'rb') as handle:
        core = MSConverter(FileContainer(FileStorage(handle)), {'mass': 230.079907196})
    assert core.typ == 'MS'
    assert core.technique is SPECTRUM_TECHNIQUES['MS']


def test_ms_composer_resolves_the_ms_descriptor_not_the_unknown_one():
    """What the descriptor buys: _technique() stops guessing from non_nmr."""
    core = JcampMSConverter(JcampBaseConverter('./tests/fixtures/source/ms/svs813f1.jdx'))
    composer = MSComposer(core)
    assert composer._technique().key == 'MS'
    assert composer._technique().x_reversed is False


def test_jcamp_ms_converter_no_longer_copies_non_nmr():
    """The property it copied is gone; the descriptor replaced it."""
    core = JcampMSConverter(JcampBaseConverter('./tests/fixtures/source/ms/svs813f1.jdx'))
    assert not hasattr(core, 'non_nmr')
    assert not hasattr(core.technique, 'non_nmr')
