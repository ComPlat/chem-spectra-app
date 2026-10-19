"""A refused file says why, in JSON — one status and one shape.

Before this, every refusal reached the user as Flask's HTML error page: 400,
403, 404 or, on most paths, an unhandled 500. That matters because of what the
ELN does with the response. `Jcamp::Create.spectrum` falls back from the
`x-extra-info-json` header to `JSON.parse(rsp.body)` and raises
`json_rsp['error']`; it only says *"Chemspectra response missing metadata
header"* when the body is **not** JSON. So every refusal, whatever its cause,
surfaced as that one sentence.

Measured on master before the change:

    /zip_jcamp        500  HTML     (fell through and returned None)
    /zip_image        500  HTML
    /jcamp            500  HTML
    /image            500  HTML
    /nmrium           500  HTML
    /zip_jcamp_n_img  403  HTML
    /combine_images   400  HTML

Nothing on the ELN side branches on the status code — every check is
`code == 200` — so one status is enough and the body carries the meaning.
"""

import io
import json

import pytest

from chem_spectra import create_app
from chem_spectra.controller.helper.refusal import REFUSED

# every route that takes a single uploaded spectrum
UPLOAD_ROUTES = [
    '/zip_jcamp_n_img',
    '/zip_jcamp',
    '/zip_image',
    '/jcamp',
    '/image',
    '/nmrium',
    '/api/v1/chemspectra/file/convert',
]


@pytest.fixture
def client():
    return create_app({'IP_WHITE_LIST': ''}).test_client()


def _post(client, route, data):
    return client.post(route, data=data, content_type='multipart/form-data')


def _refusal(response):
    """The parsed body, insisting it is JSON — which is the whole point."""
    body = response.get_data(as_text=True)
    assert body.lstrip().startswith('{'), (
        'refusal was not JSON, so the ELN will say "missing metadata header" '
        'instead of the reason: ' + body[:120])
    return json.loads(body)


# - - - nothing uploaded - - -

@pytest.mark.parametrize('route', UPLOAD_ROUTES)
def test_a_missing_file_is_refused_in_json(client, route):
    """`request.files['file']` used to raise BadRequestKeyError, which Flask
    renders as an HTML 400."""
    response = _post(client, route, {})
    assert response.status_code == REFUSED
    assert _refusal(response)['error']


@pytest.mark.parametrize('route', UPLOAD_ROUTES)
def test_an_empty_upload_is_refused_in_json(client, route):
    """An upload with no filename. `FileContainer` defines no `__bool__`, so
    the `if file:` guards were always true and this went on to fail deeper in,
    as a 500."""
    response = _post(client, route, {'file': (io.BytesIO(b''), '')})
    assert response.status_code == REFUSED
    assert 'no file' in _refusal(response)['error'].lower()


# - - - uploaded, but not convertible - - -

@pytest.mark.parametrize('route', [
    '/jcamp', '/image', '/zip_jcamp', '/zip_image', '/zip_jcamp_n_img',
    '/api/v1/chemspectra/file/convert',
])
def test_an_unconvertible_file_is_refused_in_json(client, route):
    """`/jcamp` and `/image` previously reached `send_file(False)` and raised
    `AttributeError: 'bool' object has no attribute 'read'`."""
    data = {'file': (io.BytesIO(b'not a spectrum'), 'junk.jdx')}
    response = _post(client, route, data)
    assert response.status_code == REFUSED
    assert _refusal(response)['error']


def test_the_reason_names_the_upload_not_the_code(client):
    """The message is shown to a user through the ELN, so it should describe
    the file rather than where the failure was noticed."""
    response = _post(client, '/jcamp',
                     {'file': (io.BytesIO(b'not a spectrum'), 'junk.jdx')})
    error = _refusal(response)['error']
    assert 'this file' in error.lower()
    for leak in ('Traceback', 'None', 'abort(', 'chem_spectra/'):
        assert leak not in error


# - - - combine_images - - -

def test_combine_images_refuses_with_no_files(client):
    response = _post(client, '/combine_images', {})
    assert response.status_code == REFUSED
    assert _refusal(response)['error']


# - - - the contract these all serve - - -

def test_every_refusal_uses_one_status_and_one_shape(client):
    """One convention, not several. 422 is what `UnconvertibleSpectrum`
    already returns, so refusals do not split into two kinds."""
    seen = set()
    for route in UPLOAD_ROUTES + ['/combine_images']:
        response = _post(client, route, {})
        seen.add(response.status_code)
        assert set(_refusal(response)) == {'error'}, (
            'refusal body should carry exactly one key, "error"')
    assert seen == {REFUSED}, 'refusals disagree on the status: %s' % seen
