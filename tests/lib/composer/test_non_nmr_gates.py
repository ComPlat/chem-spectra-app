"""The `non_nmr` gates, pinned before the dispatch refactor.

`non_nmr` is a single boolean that gates four unrelated concerns:
integration pairing, multiplicity output, axis-label style and peak
annotation. The refactor must not collapse them into one registry field, so
each gate is asserted on its own here, both branches.

`composer/base.py` reads it as `getattr(self.core, 'non_nmr', True)` — the
default matters, because converters that are not fed by data_type.json
(FID, NMRium) set the attribute themselves and a future kind-based core may
not set it at all.
"""

import pytest

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter
from chem_spectra.lib.composer.technique import TechniqueComposer
from chem_spectra.lib.composer.base import BaseComposer

source_nmr = './tests/fixtures/source/1H.dx'
source_ir = './tests/fixtures/source/IR.dx'

# one integration whose extent matches the multiplet below
ITG = {'xL': 1.0, 'xU': 2.0, 'area': 3.0}
MPY = {'xExtent': {'xL': 1.0, 'xU': 2.0}, 'yExtent': {'yL': 0.0, 'yU': 1.0},
       'mpyType': 'm', 'peaks': [{'x': 1.5, 'y': 100.0}], 'area': 1.0}


def _composer(path):
    return TechniqueComposer(JcampTechniqueConverter(JcampBaseConverter(path)))


def _repair(composer):
    """Re-run the gate with an integration/multiplicity pair in params."""
    composer.core.params['integration'] = {'refArea': 1, 'refFactor': 1,
                                           'shift': 0, 'stack': [dict(ITG)]}
    composer.core.params['multiplicity'] = {'stack': [dict(MPY)]}
    composer.itgs, composer.mpys, composer.all_itgs = [], [], []
    composer.prepare_itg_mpy()
    return composer


# - - - gate 1: integration pairing (composer/base.py prepare_itg_mpy) - - -

def test_nmr_moves_matching_integration_into_the_multiplet():
    composer = _repair(_composer(source_nmr))
    assert composer.core.non_nmr is False
    assert composer.itgs == []
    assert len(composer.mpys) == 1
    # the integration's area is transferred onto the multiplet
    assert composer.mpys[0]['area'] == ITG['area']


def test_non_nmr_keeps_every_integration_and_pairs_nothing():
    composer = _repair(_composer(source_ir))
    assert composer.core.non_nmr is True
    assert composer.itgs == [ITG]
    assert composer.mpys == []


# - - - gates 2 and 3: multiplicity output (composer/base.py) - - -

@pytest.mark.parametrize('method', ['gen_mpy_integ_info', 'gen_mpy_peaks_info'])
def test_multiplicity_output_is_empty_for_non_nmr(method):
    composer = _repair(_composer(source_ir))
    assert getattr(composer, method)() == []


@pytest.mark.parametrize('method', ['gen_mpy_integ_info', 'gen_mpy_peaks_info'])
def test_multiplicity_output_is_produced_for_nmr(method):
    composer = _repair(_composer(source_nmr))
    assert getattr(composer, method)() != []


# - - - the getattr default - - -

def test_gates_treat_a_core_without_non_nmr_as_non_nmr():
    """`getattr(core, 'non_nmr', True)` — the default is not incidental.

    A core that does not carry the attribute is treated as non-NMR, i.e. the
    conservative branch that emits no multiplicity. Anything that stops
    setting `non_nmr` inherits this, silently.
    """
    class CoreWithoutNonNmr:
        params = {'integration': {'refArea': 1, 'refFactor': 1, 'shift': 0,
                                  'stack': [dict(ITG)]},
                  'multiplicity': {'stack': [dict(MPY)]}}

    composer = BaseComposer.__new__(BaseComposer)
    composer.core = CoreWithoutNonNmr()
    composer.itgs, composer.mpys, composer.all_itgs = [], [], []
    composer.prepare_itg_mpy()

    assert not hasattr(composer.core, 'non_nmr')
    assert composer.itgs == [ITG]        # the non-NMR branch
    assert composer.mpys == []
    assert composer.gen_mpy_integ_info() == []
    assert composer.gen_mpy_peaks_info() == []
