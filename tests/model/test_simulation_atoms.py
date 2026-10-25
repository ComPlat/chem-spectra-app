"""`##$CSSIMULATIONATOMS`: how many atoms of the simulated nucleus the molecule
has, written beside the simulation (#339).

NMRShiftDB predicts one shift per atom it can place and leaves some out, e.g.
exchangeable NH/OH protons. The editor caps the signals a user may send for
prediction at the simulation's length; with this record it can use the
molecule's real atom count instead. The simulation itself is mocked here: the
record's content does not depend on NMRShiftDB.
"""
import io
from unittest import mock

import pytest

from fixtures.mock_predict_nmr import request_predict_nmr

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.model.molecule import MoleculeModel

ANILINE = request_predict_nmr['molfile']  # C6H7N, hydrogens implicit
SIMULATE = 'chem_spectra.model.concern.property.InferModel.simulate_nmr'


def _molecule(layout, decorate):
    return MoleculeModel(FileContainer().from_str(ANILINE), layout,
                         decorate=decorate)


@pytest.mark.parametrize('decorate', [False, True])
def test_hydrogens_are_counted_implicit_or_explicit(decorate):
    # 5 ring CH and the 2 NH protons, whether or not the molecule was built
    # with explicit hydrogens (simulation builds 1H with them)
    assert _molecule('1H', decorate).count_atoms('H') == 7


def test_carbons_are_counted():
    assert _molecule('13C', True).count_atoms('C') == 6


def _compose(client, filename, simulation):
    with open('./tests/fixtures/source/' + filename, 'rb') as f, \
            mock.patch(SIMULATE, return_value=simulation):
        response = client.post('/jcamp', data={
            'file': (f, filename),
            'molfile': (io.BytesIO(ANILINE.encode()), 'aniline.mol'),
            'simulatenmr': 'true',
        }, content_type='multipart/form-data')
    assert response.status_code == 200
    return response.data.decode()


def _record(jcamp):
    return [line for line in jcamp.splitlines()
            if line.startswith('##$CSSIMULATIONATOMS')]


@pytest.mark.parametrize('filename, simulation, atoms', [
    # NMRShiftDB leaving the two NH protons out: five shifts, seven hydrogens
    ('1H.dx', [6.6, 6.7, 6.7, 7.1, 7.1], 7),
    ('13C-CPD.dx', [115.1, 115.1, 118.5, 129.3, 129.3, 146.4], 6),
])
def test_the_record_is_written_with_the_simulation(
        client, filename, simulation, atoms):
    jcamp = _compose(client, filename, simulation)
    assert _record(jcamp) == ['##$CSSIMULATIONATOMS={}'.format(atoms)]
    # after the simulated shifts, so it does not end up inside that table
    peaks_at = jcamp.index('##$CSSIMULATIONPEAKS=')
    assert jcamp.index('##$CSSIMULATIONATOMS=') > peaks_at


def test_no_record_without_a_simulation(client):
    assert _record(_compose(client, '1H.dx', [])) == []


def test_the_record_survives_a_recompose(client):
    composed = _compose(client, '1H.dx', [6.6, 6.7, 6.7, 7.1, 7.1])
    response = client.post('/jcamp', data={
        'file': (io.BytesIO(composed.encode()), 'composed.jdx'),
    }, content_type='multipart/form-data')
    assert response.status_code == 200
    assert _record(response.data.decode()) == ['##$CSSIMULATIONATOMS=7']
