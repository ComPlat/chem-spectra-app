"""Issue #244 — two valid BagIt archives in one upload, only one processed.

Reported 2025-08-06 with a video and never investigated. Reproduced here:
a zip whose top level is two BagIt directories returned **3 curves**, the
contents of whichever `os.walk` reached first. The other archive's 2 curves
were dropped with no message.

The cause is one line. `find_dir` (`model/transformer.py`) returns the first
directory containing `bagit.txt` and stops, and `to_composer` hands that single
directory to `BagItBaseConverter`.

Processing both was the wrong fix. The reporter said so in the issue: *"two
bagits are normally different Datasets"*, and *"ELN needs to send feedback /
warning to user, to ask them to create additional Datasets"*. Merging two
datasets into one attachment group would be incorrect, not merely surprising.
So this refuses, and says how many it found.
"""

import io
import zipfile

import pytest

from chem_spectra import create_app

BAGIT_A = './tests/fixtures/source/bagit/cv/File053_BagIt.zip'   # 3 curves
BAGIT_B = './tests/fixtures/source/bagit/aif/aif.zip'            # 2 curves


@pytest.fixture
def client():
    return create_app({'IP_WHITE_LIST': ''}).test_client()


def _nested(*sources):
    """One zip whose top level is a directory per BagIt archive."""
    out = io.BytesIO()
    with zipfile.ZipFile(out, 'w') as archive:
        for index, source in enumerate(sources):
            with zipfile.ZipFile(source) as inner:
                for name in inner.namelist():
                    archive.writestr('bag_%d/%s' % (index, name), inner.read(name))
    out.seek(0)
    return out


def _post(client, payload):
    return client.post(
        '/zip_jcamp_n_img',
        data={'file': (payload, 'upload.zip')},
        content_type='multipart/form-data',
    )


def _curves(source):
    with zipfile.ZipFile(source) as archive:
        return len([n for n in archive.namelist()
                    if n.startswith('data/') and n.lower().endswith('.jdx')])


def test_two_bagits_are_refused_rather_than_silently_halved(client):
    response = _post(client, _nested(BAGIT_A, BAGIT_B))
    assert response.status_code == 422
    error = response.get_json()['error']
    assert '2 BagIt archives' in error
    assert 'one at a time' in error


def test_the_refusal_says_how_many_were_found(client):
    """Three, so the count is not a coincidence of the two-archive case."""
    response = _post(client, _nested(BAGIT_A, BAGIT_B, BAGIT_A))
    assert response.status_code == 422
    assert '3 BagIt archives' in response.get_json()['error']


def test_one_bagit_is_unaffected(client):
    """The guard must not cost the ordinary case. This is what used to come
    back for the two-archive upload as well, which is how the loss hid."""
    with open(BAGIT_A, 'rb') as handle:
        response = _post(client, io.BytesIO(handle.read()))
    assert response.status_code == 200
    with zipfile.ZipFile(io.BytesIO(response.get_data())) as archive:
        produced = [n for n in archive.namelist() if n.lower().endswith('.jdx')]
    assert len(produced) == _curves(BAGIT_A) == 3


def test_a_nested_single_bagit_still_works(client):
    """One archive inside a wrapping directory — the shape that made the two
    case look like a normal upload."""
    response = _post(client, _nested(BAGIT_A))
    assert response.status_code == 200


def test_find_dirs_returns_every_match_not_the_first():
    """The unit beneath it, so a regression names the cause rather than the
    symptom."""
    import os
    import tempfile
    from chem_spectra.model.transformer import find_dirs, find_dir

    with tempfile.TemporaryDirectory() as root:
        for name in ('alpha', 'beta'):
            os.makedirs(os.path.join(root, name))
            open(os.path.join(root, name, 'bagit.txt'), 'w').close()
        assert len(find_dirs(root, 'bagit.txt')) == 2
        assert find_dirs(root, 'bagit.txt') == sorted(find_dirs(root, 'bagit.txt'))
        # find_dir keeps its first-match behaviour, which Bruker's `fid` needs
        assert find_dir(root, 'bagit.txt') in find_dirs(root, 'bagit.txt')
