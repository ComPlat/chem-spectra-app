"""Request-size and archive limits, and how a request over them is refused.

Every refusal is the usual 422 with `{"error": ...}` naming the sizes, and on
the prediction routes the usual HTTP 200 with the reason in `outline`.
"""
import io
import zipfile

import pytest
from werkzeug.datastructures import FileStorage

from chem_spectra import DEFAULT_MAX_CONTENT_LENGTH, create_app
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.lib.converter.share import UnconvertibleSpectrum


SOURCE = './tests/fixtures/source/'
MB = 1024 * 1024


def _read(name):
    with open(SOURCE + name, 'rb') as f:
        return f.read()


def _zip(members):
    """An in-memory zip of {name: bytes}."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as z:
        for name, data in members.items():
            z.writestr(name, data)
    return buf.getvalue()


def _post(client, url, **files):
    data = {k: (io.BytesIO(v), name) for k, (name, v) in files.items()}
    return client.post(url, data=data, content_type='multipart/form-data')


# -- the request body --------------------------------------------------------

def test_the_default_request_limit_is_200_mb(app):
    assert app.config['MAX_CONTENT_LENGTH'] == DEFAULT_MAX_CONTENT_LENGTH
    assert DEFAULT_MAX_CONTENT_LENGTH == 200 * MB


def test_the_request_limit_can_be_set_from_the_environment(monkeypatch):
    monkeypatch.setenv('FLASK_MAX_CONTENT_LENGTH', str(5 * MB))
    assert create_app().config['MAX_CONTENT_LENGTH'] == 5 * MB


def test_an_upload_over_the_limit_is_refused_with_the_sizes(app, client):
    app.config['MAX_CONTENT_LENGTH'] = MB
    response = _post(client, '/zip_jcamp_n_img',
                     file=('big.jdx', b'#' * (2 * MB)))
    assert response.status_code == 422
    assert response.mimetype == 'application/json'
    error = response.json['error']
    assert error.startswith('the upload is 2.0 MB')
    assert 'the limit is 1.0 MB' in error


def test_an_upload_under_the_limit_is_processed(app, client):
    app.config['MAX_CONTENT_LENGTH'] = MB
    response = _post(client, '/zip_jcamp_n_img', file=('IR.dx', _read('IR.dx')))
    assert response.status_code == 200


def test_a_form_field_over_its_limit_is_refused(app, client):
    app.config['MAX_FORM_MEMORY_SIZE'] = 1000
    response = client.post(
        '/zip_jcamp_n_img',
        data={'file': (io.BytesIO(_read('IR.dx')), 'IR.dx'),
              'peaks_str': 'x' * 5000},
        content_type='multipart/form-data')
    assert response.status_code == 422
    assert 'MAX_FORM_MEMORY_SIZE' in response.json['error']


def test_a_prediction_over_the_limit_keeps_the_outline(app, client):
    app.config['MAX_CONTENT_LENGTH'] = MB
    response = _post(client, '/api/v1/chemspectra/predict/infrared',
                     spectrum=('big.jdx', b'#' * (2 * MB)),
                     molfile=('m.mol', _read('molfile/svs813f1_B.mol')))
    assert response.status_code == 200
    outline = response.json['outline']
    assert outline['code'] == 422
    assert outline['text'].startswith('the upload is 2.0 MB')


# -- archives ----------------------------------------------------------------

def test_an_archive_over_the_size_limit_is_refused_with_the_sizes(app, client):
    app.config['MAX_ZIP_SIZE'] = 1
    archive = _zip({'fid': b'\0' * (2 * MB)})  # compresses to a few kB
    for url in ('/zip_jcamp_n_img', '/api/v1/chemspectra/file/convert'):
        response = _post(client, url, file=('a.zip', archive))
        assert response.status_code == 422, url
        assert response.json['error'] == \
            'this archive expands to 2.0 MB; the limit is 1 MB.'


def test_an_archive_over_the_member_limit_is_refused(app, client):
    app.config['MAX_ZIP_MEMBERS'] = 3
    archive = _zip({'f{}'.format(i): b'x' for i in range(5)})
    response = _post(client, '/zip_jcamp_n_img', file=('a.zip', archive))
    assert response.status_code == 422
    assert response.json['error'] == \
        'this archive holds 5 files; the limit is 3.'


def test_a_prediction_with_an_archive_over_the_limit_keeps_the_outline(
        app, client):
    app.config['MAX_ZIP_SIZE'] = 1
    response = _post(client, '/api/v1/chemspectra/predict/infrared',
                     spectrum=('a.zip', _zip({'fid': b'\0' * (2 * MB)})),
                     molfile=('m.mol', _read('molfile/svs813f1_B.mol')))
    assert response.status_code == 200
    assert response.json['outline'] == {
        'code': 422,
        'text': 'this archive expands to 2.0 MB; the limit is 1 MB.',
    }


def test_the_largest_fixture_archive_is_within_the_defaults():
    # lcms.zip unpacks to ~120 MB across 7 files; the defaults are 200 MB and
    # 1000 files. Outside an app context the defaults apply.
    with open(SOURCE + 'lcms/lcms.zip', 'rb') as f:
        container = FileContainer(FileStorage(f, filename='lcms.zip'))
    assert container.bcore


def test_the_limits_apply_outside_a_request_too():
    archive = _zip({'f{}'.format(i): b'x' for i in range(1001)})
    with pytest.raises(UnconvertibleSpectrum, match='holds 1001 files'):
        FileContainer(FileStorage(io.BytesIO(archive), filename='a.zip'))
