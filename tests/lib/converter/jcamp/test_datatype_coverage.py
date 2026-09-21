"""Every datatype chemotion-converter-app can emit is accounted for here.

The two lists live in different repositories and nothing keeps them in
step. That drift is what produced PR #291: twelve datatypes were unmapped
in `data_type.json`, five of them auxiliary and correctly so, the rest not
-- and files of those types were silently misclassified.

A direct import is impossible across repos, so this compares against a
checked-in snapshot, `tests/fixtures/converter_app_data_types.txt`. That
only catches drift when the snapshot is refreshed, so it is a backstop and
not a substitute for being told when `DATA_TYPES` changes upstream.

To refresh the snapshot, from a checkout of chemotion-converter-app::

    git -C <converter-app> show origin/master:converter_app/options.py \
      | sed -n '/^DATA_TYPES = (/,/^)/p'

Take the values from `origin/master`, not from whatever branch happens to
be checked out -- when this was written the local checkout was on a feature
branch five commits behind and missing two datatypes.
"""

import json
import os
import re

import pytest

from chem_spectra.lib.converter.jcamp import base as base_module

SNAPSHOT = os.path.join(
    os.path.dirname(__file__), '..', '..', '..', 'fixtures',
    'converter_app_data_types.txt',
)

# Datatypes that must NOT be mapped, each with the reason. A datatype in
# neither this set nor data_type.json is what this module exists to catch.
INTENTIONALLY_UNMAPPED = {
    # --- auxiliary blocks -------------------------------------------------
    # These sit alongside a primary block in the same file. Block selection
    # takes the first *recognised* datatype, so mapping one of these would
    # make every affected file read the auxiliary block instead of its
    # spectrum. tests/lib/converter/jcamp/test_jcamp_datatype_classification
    # .py::test_auxiliary_blocks_stay_unmapped enforces this directly.
    'NMR FID': 'auxiliary block; mapping it would read the FID, not the spectrum',
    'NMR PEAK TABLE': 'auxiliary block',
    'NMP PEAK ASSIGNMENTS': 'auxiliary block (upstream misspelling of NMR)',
    'INFRARED PEAK TABLE': 'auxiliary block',
    'INFRARED INTERFEROGRAM': 'auxiliary block; the raw interferogram, not the spectrum',

    # --- techniques deferred by decision ----------------------------------
    # Since #291 these take the generic curve path with a warning rather
    # than crashing, so adding them is not urgent. Each needs axis
    # conventions, a threshold and a plotting decision -- domain calls.
    'SQUID': 'no axis conventions decided; magnetometry',
    'TENSIOMETRY': 'no axis conventions decided',
    'LINEAR SWEEP VOLTAMMETRY': 'must not silently borrow the CV path, which assumes a cyclic sweep',
    'SINGLE CRYSTAL X-RAY DIFFRACTION': 'different measurement from powder XRD, different plot conventions',
    'INFRARED TRANSFERED SPECTRUM': 'unclear whether it is a primary spectrum',
}


def _snapshot_datatypes():
    with open(SNAPSHOT, encoding='utf-8') as handle:
        return [line.strip() for line in handle
                if line.strip() and not line.startswith('#')]


def _mapped_values():
    path = os.path.join(os.path.dirname(base_module.__file__), 'data_type.json')
    with open(path, encoding='utf-8') as handle:
        mapping = json.load(handle)['datatypes']
    return {value.upper() for values in mapping.values() for value in values}


def test_snapshot_is_readable_and_attributed():
    """A snapshot with no provenance cannot be refreshed responsibly."""
    with open(SNAPSHOT, encoding='utf-8') as handle:
        header = handle.read().split('\n\n')[0]
    assert re.search(r'# Commit : [0-9a-f]{12}', header)
    assert '# Taken  :' in header
    assert len(_snapshot_datatypes()) > 20


@pytest.mark.parametrize('datatype', _snapshot_datatypes())
def test_every_upstream_datatype_is_mapped_or_deliberately_not(datatype):
    """The whole point: a new upstream datatype fails here by name.

    If this fails, the datatype is new since the snapshot was taken. Either
    add it to data_type.json with a SPECTRUM_TECHNIQUES entry, or add it to
    INTENTIONALLY_UNMAPPED above with a one-line reason. Do not delete the
    row from the snapshot.
    """
    mapped = datatype.upper() in _mapped_values()
    deliberate = datatype in INTENTIONALLY_UNMAPPED
    assert mapped or deliberate, (
        '%r is emitted by chemotion-converter-app but is neither mapped in '
        'data_type.json nor listed in INTENTIONALLY_UNMAPPED. Map it, or '
        'record why it should not be.' % datatype
    )


def test_a_datatype_is_not_both_mapped_and_excluded():
    """Contradiction between the mapping and the exclusion list."""
    both = {d for d in INTENTIONALLY_UNMAPPED if d.upper() in _mapped_values()}
    assert both == set(), (
        'listed as intentionally unmapped but present in data_type.json: %s'
        % sorted(both)
    )


def test_exclusions_all_still_exist_upstream():
    """A stale exclusion hides the fact that upstream dropped a datatype."""
    stale = set(INTENTIONALLY_UNMAPPED) - set(_snapshot_datatypes())
    assert stale == set(), (
        'listed as intentionally unmapped but absent from the upstream '
        'snapshot -- upstream may have removed them: %s' % sorted(stale)
    )


def test_every_exclusion_carries_a_reason():
    for datatype, reason in INTENTIONALLY_UNMAPPED.items():
        assert reason and len(reason) > 10, datatype
