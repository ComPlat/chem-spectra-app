"""The spectrum-technique registry must agree with data_type.json and with the
behaviour the pre-refactor code exhibits.

The parity test is the check that stops these two drifting apart the way
data_type.json drifted twelve entries behind chemotion-converter-app.
"""

import json
import os

import pytest

from chem_spectra.lib.converter.jcamp import base as base_module
from chem_spectra.lib.converter.jcamp.techniques import (
    SPECTRUM_TECHNIQUES, UNKNOWN_TECHNIQUE, SpectrumTechnique, technique_for,
)
from tests.lib.composer.test_technique_rendering_matrix import EXPECTED


def _mapping_keys():
    path = os.path.join(os.path.dirname(base_module.__file__), 'data_type.json')
    with open(path) as handle:
        return set(json.load(handle)['datatypes'])


def test_every_mapped_datatype_has_a_technique():
    assert _mapping_keys() - set(SPECTRUM_TECHNIQUES) == set()


def test_every_technique_corresponds_to_a_mapped_datatype():
    assert set(SPECTRUM_TECHNIQUES) - _mapping_keys() == set()


def test_each_technique_knows_its_own_key():
    for key, technique in SPECTRUM_TECHNIQUES.items():
        assert technique.key == key


@pytest.mark.parametrize('key', sorted(EXPECTED))
def test_technique_reproduces_the_measured_behaviour_matrix(key):
    """The registry must agree with what the current code actually draws.

    EXPECTED is imported from the rendering matrix rather than restated, so
    the two cannot disagree silently.
    """
    threshold, orientation, x_style, y_style = EXPECTED[key]
    technique = SPECTRUM_TECHNIQUES[key]

    assert technique.threshold == threshold
    assert technique.x_reversed is (orientation == 'reversed')
    assert technique.x_axis == x_style
    assert technique.y_axis == y_style


def test_unknown_technique_matches_the_generic_curve_path():
    assert technique_for('NEUTRON SCATTERING') is UNKNOWN_TECHNIQUE
    assert technique_for('') is UNKNOWN_TECHNIQUE
    assert UNKNOWN_TECHNIQUE.threshold == 0.5
    assert UNKNOWN_TECHNIQUE.x_reversed is True
    assert UNKNOWN_TECHNIQUE.x_axis == 'generic'
    assert UNKNOWN_TECHNIQUE.y_axis == 'generic'


def test_only_nmr_carries_the_nmr_gates():
    """All four are True only for NMR today.

    They are separate fields so that a future technique can enable one
    without the others; this test records that today they move together.
    """
    for key, technique in SPECTRUM_TECHNIQUES.items():
        gates = (technique.nmr_integration, technique.multiplicity,
                 technique.peak_annotation, technique.nmr_headers)
        assert all(gates) if key == 'NMR' else not any(gates), key


def test_em_wave_grouping_matches_the_old_predicate():
    """`is_em_wave` was `typ in ['INFRARED', 'RAMAN', 'UVVIS']`."""
    grouped = {k for k, technique in SPECTRUM_TECHNIQUES.items() if technique.em_wave}
    assert grouped == {'INFRARED', 'RAMAN', 'UVVIS'}


def test_cyclic_voltammetry_flag_is_cyclic_voltammetry_only():
    """Formerly cv_scaling, which nothing consumed. All 14 CV sites read it now."""
    flagged = {
        k for k, technique in SPECTRUM_TECHNIQUES.items()
        if technique.cyclic_voltammetry
    }
    assert flagged == {'CYCLIC VOLTAMMETRY'}


def test_techniques_are_immutable():
    with pytest.raises(Exception):
        SPECTRUM_TECHNIQUES['NMR'].threshold = 1.0
