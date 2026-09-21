"""A JCAMP nmrglue cannot extract a data array from.

`make_ni_data_ys` dereferenced `base.data.shape` unguarded. nmrglue returns
None when it parses no data array, so such a file raised AttributeError
straight out of the request -- a 500 for what is really just an unusable
upload.

It now raises UnparsableJcampData, which jcamp2cvp turns into the same
falsy result any other unconvertible file produces, so the controllers'
existing `abort(...)` handles it. No new contract was invented: this is the
status those endpoints already return when conversion yields nothing.
"""

import io

import pytest

from chem_spectra import create_app
from chem_spectra.lib.converter.jcamp.data_parse import UnparsableJcampData

# 129 data-block markers, but nmrglue yields data=None
source_unparsable = './tests/fixtures/source/MS.dx'


@pytest.fixture
def client():
    return create_app({'TESTING': True}).test_client()


def _post(client, endpoint):
    with open(source_unparsable, 'rb') as handle:
        return client.post(
            endpoint,
            data={'file': (io.BytesIO(handle.read()), 'unparsable.jdx')},
            content_type='multipart/form-data',
        )


@pytest.mark.parametrize('endpoint,expected', [
    # each endpoint's own existing "could not convert" status, unchanged
    ('/api/v1/chemspectra/file/convert', 400),
    ('/zip_jcamp_n_img', 403),
])
def test_unparsable_jcamp_is_rejected_not_a_500(client, endpoint, expected):
    """Asserted at the endpoint. A converter-level test would pass while
    the controller still returned 500 -- that has happened twice."""
    assert _post(client, endpoint).status_code == expected


def test_the_parse_failure_is_named(client):
    """The exception says what went wrong, rather than AttributeError on
    'NoneType' from three frames down."""
    from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
    from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter

    base = JcampBaseConverter(source_unparsable)
    assert base.data is None
    with pytest.raises(UnparsableJcampData):
        JcampNIConverter(base)


def test_a_parsable_file_is_unaffected(client):
    """The guard must not catch working files."""
    with open('./tests/fixtures/source/IR.dx', 'rb') as handle:
        response = client.post(
            '/api/v1/chemspectra/file/convert',
            data={'file': (io.BytesIO(handle.read()), 'IR.dx')},
            content_type='multipart/form-data',
        )
    assert response.status_code == 200
