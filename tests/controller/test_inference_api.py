import pytest
import io
import json
from unittest import mock
from fixtures.mock_predict_nmr import RequestPredictNmr, ResponsePredictNmr
from fixtures.mock_predict_ir import ResponsePredictIr


target_dir = './tests/fixtures/'
source_dir = 'source/'


@mock.patch('requests.post', mock.Mock(return_value=ResponsePredictNmr()))
def test_api_chemspectra_predict_by_peaks_json(client):
    response = client.post(
        '/api/v1/chemspectra/predict/by_peaks_json',
        content_type='application/json',
        data=RequestPredictNmr().json()
    )

    assert response.status_code == 200
    assert response.json['status'] == True
    assert response.mimetype == 'application/json'


@mock.patch('requests.post', mock.Mock(return_value=ResponsePredictNmr()))
def test_api_chemspectra_predict_by_peaks_form(client):
    with open(target_dir + source_dir + '/molfile/svs813f1_B.mol', 'rb') as f:
        file_content = f.read()
    params = RequestPredictNmr().json()
    params = json.loads(params)
    data = dict(
        molfile=(io.BytesIO(file_content), 'svs813f1_B.mol'),
        layout=params['layout'],
        peaks=json.dumps(params['peaks']),
        shift=json.dumps(params['shift']),
    )
    response = client.post(
        '/api/v1/chemspectra/predict/by_peaks_form',
        content_type='multipart/form-data',
        data=data
    )
    assert response.status_code == 200
    assert response.json['status'] == True
    assert response.mimetype == 'application/json'


@mock.patch('requests.post', mock.Mock(return_value=ResponsePredictIr()))
def test_api_chemspectra_predict_infrared(client):
    with open(target_dir + source_dir + '/IR.dx', 'rb') as f:
        spectrum_content = f.read()
    with open(target_dir + source_dir + '/molfile/svs813f1_B.mol', 'rb') as f:
        molfile_content = f.read()
    data = {}
    data['molfile'] = (io.BytesIO(molfile_content), 'svs813f1_B.mol')
    data['spectrum'] = (io.BytesIO(spectrum_content), 'IR.dx')

    response = client.post(
        '/api/v1/chemspectra/predict/infrared',
        content_type='multipart/form-data',
        data=data
    )
    assert response.status_code == 200
    assert json.loads(response.json)['status'] == True
    assert response.mimetype == 'application/json'
