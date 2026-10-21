"""Every prediction failure is HTTP 200 with a reason in `outline`.

The ELN keeps a prediction response only when the status is 200, and the
editor's forecast panel reads only `outline.code`. Before this, a missing
upload was Flask's HTML 400, an unsupported layout or an unset predictor URL a
500, and a predictor that never answered held the worker until gunicorn
killed it. None of that reached the user.
"""
import io
import json
from unittest import mock

import pytest
import requests

from fixtures.mock_predict_nmr import request_predict_nmr, response_predict_nmr


SOURCE = './tests/fixtures/source/'
NMR_FORM = '/api/v1/chemspectra/predict/nmr_peaks_form'
NMR_JSON = '/api/v1/chemspectra/predict/nmr_peaks_json'
IR = '/api/v1/chemspectra/predict/infrared'
MS = '/api/v1/chemspectra/predict/ms'


class Answer:
    def __init__(self, status_code=200, body=None, text=None):
        self.status_code = status_code
        self._body = body
        self._text = text

    def json(self):
        if self._text is not None:
            raise requests.exceptions.JSONDecodeError('not JSON', self._text, 0)
        return self._body


@pytest.fixture
def predictors(app):
    app.config['URL_NSHIFTDB'] = 'http://nmrshiftdb.test/quickcheck'
    app.config['URL_DEEPIR'] = 'http://deepir.test/infer_ir'


def _read(name):
    with open(SOURCE + name, 'rb') as f:
        return f.read()


def _molfile():
    return (io.BytesIO(_read('molfile/svs813f1_B.mol')), 'svs813f1_B.mol')


def _nmr_form(**overrides):
    data = dict(
        molfile=_molfile(),
        layout=request_predict_nmr['layout'],
        peaks=json.dumps(request_predict_nmr['peaks']),
        shift=json.dumps(request_predict_nmr['shift']),
    )
    data.update(overrides)
    return {k: v for k, v in data.items() if v is not None}


def _post(client, url, data):
    return client.post(url, content_type='multipart/form-data', data=data)


def _outline(response):
    assert response.status_code == 200
    assert response.mimetype == 'application/json'
    return response.json['outline']


# -- requests the server cannot use: 422, no predictor asked ------------------

@pytest.mark.parametrize('overrides, reason', [
    ({'molfile': None}, 'molfile'),
    ({'molfile': (io.BytesIO(b'not a molfile'), 'x.mol')}, 'molfile'),
    ({'peaks': '{}'}, 'peaks'),
    ({'peaks': 'not json'}, 'JSON'),
    ({'layout': 'IR'}, 'NMR prediction covers 1H, 13C, 19F, not IR'),
])
def test_nmr_form_refuses_what_it_cannot_ask_about(
        client, predictors, overrides, reason):
    with mock.patch('requests.post') as post:
        outline = _outline(_post(client, NMR_FORM, _nmr_form(**overrides)))
    assert outline['code'] == 422
    assert reason in outline['text']
    post.assert_not_called()


@pytest.mark.parametrize('body, reason', [
    ('not json', 'JSON object'),
    (json.dumps({**request_predict_nmr, 'molfile': None}), 'molfile'),
    (json.dumps({**request_predict_nmr, 'peaks': []}), 'peaks'),
])
def test_nmr_json_refuses_what_it_cannot_ask_about(
        client, predictors, body, reason):
    with mock.patch('requests.post') as post:
        response = client.post(
            NMR_JSON, content_type='application/json', data=body)
    outline = _outline(response)
    assert outline['code'] == 422
    assert reason in outline['text']
    post.assert_not_called()


def test_ir_refuses_a_request_without_a_spectrum(client, predictors):
    outline = _outline(_post(client, IR, {'molfile': _molfile()}))
    assert outline['code'] == 422
    assert 'spectrum' in outline['text']


def test_ms_refuses_a_request_without_a_spectrum(client, predictors):
    outline = _outline(_post(client, MS, {'molfile': _molfile()}))
    assert outline['code'] == 422
    assert 'spectrum' in outline['text']


# -- the predictor is not there: 503 ----------------------------------------

def test_nmr_without_a_configured_predictor(client):
    with mock.patch('requests.post') as post:
        outline = _outline(_post(client, NMR_FORM, _nmr_form()))
    assert outline == {
        'code': 503,
        'text': 'NMR prediction is not configured on this server.',
    }
    post.assert_not_called()


def test_ir_without_a_configured_predictor(client):
    # production today: bakery's config leaves URL_DEEPIR commented out
    with mock.patch('requests.post') as post:
        outline = _outline(_post(client, IR, {
            'molfile': _molfile(),
            'spectrum': (io.BytesIO(_read('IR.dx')), 'IR.dx'),
        }))
    assert outline == {
        'code': 503,
        'text': 'IR prediction is not configured on this server.',
    }
    post.assert_not_called()


def test_nmr_predictor_unreachable(client, predictors):
    with mock.patch('requests.post', side_effect=requests.ConnectionError):
        outline = _outline(_post(client, NMR_FORM, _nmr_form()))
    assert outline['code'] == 503


# -- the predictor does not answer, or answers badly ------------------------

def test_nmr_predictor_timeout(client, app, predictors):
    app.config['PREDICT_TIMEOUT'] = 7
    with mock.patch('requests.post', side_effect=requests.Timeout) as post:
        outline = _outline(_post(client, NMR_FORM, _nmr_form()))
    assert outline == {
        'code': 504, 'text': 'NMRShiftDB did not answer within 7 s.'}
    assert post.call_args.kwargs['timeout'] == 7


def test_ir_predictor_timeout(client, app, predictors):
    app.config['PREDICT_TIMEOUT'] = 7
    with mock.patch('requests.post', side_effect=requests.Timeout) as post:
        outline = _outline(_post(client, IR, {
            'molfile': _molfile(),
            'spectrum': (io.BytesIO(_read('IR.dx')), 'IR.dx'),
        }))
    assert outline['code'] == 504
    assert post.call_args.kwargs['timeout'] == 7


def test_both_predictors_get_a_timeout_by_default(client, predictors):
    with mock.patch('requests.post', return_value=Answer(
            body=response_predict_nmr)) as post:
        _post(client, NMR_FORM, _nmr_form())
    assert post.call_args.kwargs['timeout'] == 120


@pytest.mark.parametrize('answer, code', [
    (Answer(status_code=503, text='<html>down</html>'), 502),
    (Answer(body={'unexpected': True}), 502),
    (Answer(body={'result': []}), 502),
    # a body that is not JSON stays the "Peak & Element count mismatch" it
    # was before
    (Answer(text='Number of peaks exceeds number of atoms'), 400),
])
def test_nmr_predictor_bad_answers(client, predictors, answer, code):
    with mock.patch('requests.post', return_value=answer):
        outline = _outline(_post(client, NMR_FORM, _nmr_form()))
    assert outline['code'] == code


@pytest.mark.parametrize('answer', [
    Answer(status_code=500, text='<html>boom</html>'),
    Answer(text='<html>not json</html>'),
    Answer(body={'outline': {'code': 200}}),
])
def test_ir_predictor_bad_answers(client, predictors, answer):
    with mock.patch('requests.post', return_value=answer):
        outline = _outline(_post(client, IR, {
            'molfile': _molfile(),
            'spectrum': (io.BytesIO(_read('IR.dx')), 'IR.dx'),
        }))
    assert outline['code'] == 502
