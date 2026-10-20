"""Every fixture still produces what it is recorded as producing.

Each fixture is read and composed the way the endpoints do it, and its
observables -- the chosen block, x and y, factors, labels, peaks, integrals,
the MS runs and the composed JCAMP's hash -- are compared with its record
under `tests/snapshots/`. A failure names the fixture and shows which
observable moved.

When a change is meant to alter output, re-record and commit the records with
the change, so the diff of `tests/snapshots/` says what moved:

    python -m tests.tools.snapshot_baseline
"""

import json

import pytest

from tests.tools import snapshot_baseline as harness

FIXTURES = harness.fixtures()


@pytest.mark.parametrize(
    'path', FIXTURES, ids=[str(p.relative_to(harness.FIXTURES)) for p in FIXTURES])
def test_the_fixture_produces_what_it_is_recorded_as_producing(path):
    stored = harness.snapshot_path(path)
    assert stored.exists(), (
        'no record for %s; run `python -m tests.tools.snapshot_baseline`'
        % path.name)
    assert json.loads(harness.render(path)) == json.loads(
        stored.read_text(encoding='utf-8'))
