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

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.lib.converter.ms import MSConversionFailed, MSConverterUnavailable
from chem_spectra.model.inferencer import InferencerModel as InferModel
from chem_spectra.model.molecule import MoleculeModel


SOURCE = './tests/fixtures/source/'
NMR_FORM = '/api/v1/chemspectra/predict/nmr_peaks_form'
NMR_JSON = '/api/v1/chemspectra/predict/nmr_peaks_json'
IR = '/api/v1/chemspectra/predict/infrared'
MS = '/api/v1/chemspectra/predict/ms'


class Answer:
    """A predictor's answer: a JSON `body`, or a `text` that is not JSON."""

    def __init__(self, status_code=200, body=None, text=None):
        self.status_code = status_code
        self._body = body
        self.text = json.dumps(body) if text is None else text

    def json(self):
        if self._body is None:
            raise requests.exceptions.JSONDecodeError('not JSON', self.text, 0)
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


# -- review of 2026-10-22: inputs that still reached a 500 ------------------

@pytest.mark.parametrize('overrides, reason', [
    ({'peaks': '{"a": 1}'}, 'list'),
    ({'peaks': '[{"y": 1}]'}, 'numeric "x"'),
    ({'peaks': '[{"x": "7.2"}]'}, 'numeric "x"'),
    ({'peaks': 'null'}, 'No peaks'),
    ({'shift': '[1]'}, 'shift must be a JSON object'),
])
def test_nmr_form_refuses_malformed_peaks_and_shift(
        client, predictors, overrides, reason):
    with mock.patch('requests.post') as post:
        outline = _outline(_post(client, NMR_FORM, _nmr_form(**overrides)))
    assert outline['code'] == 422
    assert reason in outline['text']
    post.assert_not_called()


@pytest.mark.parametrize('shift', ['null', '{"ref": null}', '{"ref": "x"}'])
def test_nmr_form_takes_a_shift_without_a_usable_ref(client, predictors, shift):
    with mock.patch('requests.post', return_value=Answer(
            body=response_predict_nmr)) as post:
        outline = _outline(_post(client, NMR_FORM, _nmr_form(shift=shift)))
    assert outline['code'] == 200
    assert post.call_args.kwargs['json']['inputs'][0]['solvent'] == 'None'


def test_nmr_json_refuses_peaks_that_are_not_a_list(client, predictors):
    body = json.dumps({**request_predict_nmr, 'peaks': {'x': 1}})
    with mock.patch('requests.post') as post:
        response = client.post(
            NMR_JSON, content_type='application/json', data=body)
    assert _outline(response)['code'] == 422
    post.assert_not_called()


def _ir(client):
    return _post(client, IR, {
        'molfile': _molfile(),
        'spectrum': (io.BytesIO(_read('IR.dx')), 'IR.dx'),
    })


@pytest.mark.parametrize('fgs', [
    ['C=O'],
    [{'smarts': 'C=O'}],
    [{'sma': 'C(=O'}],
    'C=O',
])
def test_ir_answer_with_groups_that_cannot_be_drawn(client, predictors, fgs):
    answer = Answer(body={'outline': {'code': 200}, 'output': {
        'result': [{'type': 'ir', 'fgs': fgs}]}})
    with mock.patch('requests.post', return_value=answer):
        outline = _outline(_ir(client))
    assert outline['code'] == 502


def test_ms_refuses_a_spectrum_that_is_not_a_mass_spectrum(client, predictors):
    outline = _outline(_post(client, MS, {
        'molfile': _molfile(),
        'spectrum': (io.BytesIO(_read('1H.dx')), '1H.dx'),
    }))
    assert outline == {'code': 422, 'text': 'MS prediction needs a mass spectrum.'}


@pytest.mark.parametrize('error, code', [
    (MSConverterUnavailable('msconvert did not answer'), 503),
    (MSConversionFailed('no mzML for this upload'), 422),
])
def test_ms_conversion_failure_stays_in_the_outline(
        client, predictors, error, code):
    with mock.patch(
            'chem_spectra.controller.inference_api.TraModel.to_composer',
            side_effect=error):
        outline = _outline(_post(client, MS, {
            'molfile': _molfile(),
            'spectrum': (io.BytesIO(b'raw'), 'x.raw'),
        }))
    assert outline['code'] == code
    assert str(error) in outline['text']


def test_13c_conversion_failure_stays_in_the_outline(client, predictors):
    with mock.patch(
            'chem_spectra.controller.inference_api.TraModel.to_converter',
            side_effect=MSConverterUnavailable('msconvert did not answer')):
        outline = _outline(_post(client, NMR_FORM, _nmr_form(
            peaks='{}', spectrum=(io.BytesIO(b'x'), 'x.mzML'))))
    assert outline['code'] == 503


# -- review of 2026-10-22: answers that lost the predictor's reason ---------

def test_ir_predictor_keeps_its_own_reason(client, predictors):
    answer = Answer(status_code=400, body={
        'outline': {'code': 400, 'text': 'unsupported molecule'}})
    with mock.patch('requests.post', return_value=answer):
        outline = _outline(_ir(client))
    assert outline == {'code': 400, 'text': 'unsupported molecule'}


def test_nmr_error_body_is_quoted(client, predictors):
    answer = Answer(status_code=400, body={'error': 'bad solvent'})
    with mock.patch('requests.post', return_value=answer):
        outline = _outline(_post(client, NMR_FORM, _nmr_form()))
    assert outline['code'] == 502
    assert 'HTTP 400' in outline['text']
    assert 'bad solvent' in outline['text']


def test_ir_server_error_is_reported_as_an_answer(client, predictors):
    with mock.patch('requests.post', return_value=Answer(
            status_code=500, text='<html>boom</html>')):
        outline = _outline(_ir(client))
    assert outline['code'] == 502
    assert 'answered HTTP 500' in outline['text']
    assert 'could not be asked' not in outline['text']


def test_simulation_uses_predict_timeout(app, predictors):
    app.config['PREDICT_TIMEOUT'] = 7
    with app.app_context():
        mm = MoleculeModel(FileContainer().from_str(request_predict_nmr['molfile']), '13C')
        with mock.patch('requests.post', side_effect=requests.Timeout) as post:
            assert InferModel.simulate_nmr(mm=mm, layout='13C') == []
    assert post.call_args.kwargs['timeout'] == 7
