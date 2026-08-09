import io
import json
from unittest import mock
from fixtures.mock_predict_nmr import RequestPredictNmr, ResponsePredictNmr
from fixtures.mock_predict_ir import ResponsePredictIr

from tests.dataset_catalog import dataset_path


@mock.patch('requests.post', mock.Mock(return_value=ResponsePredictNmr()))
def test_api_chemspectra_predict_by_peaks_json(client):
    response = client.post(
        '/api/v1/chemspectra/predict/nmr_peaks_json',
        content_type='application/json',
        data=RequestPredictNmr().json()
    )
    assert response.status_code == 200
    assert response.json['outline']['code'] == 200
    assert response.mimetype == 'application/json'


@mock.patch('requests.post', mock.Mock(return_value=ResponsePredictNmr()))
def test_api_chemspectra_predict_by_peaks_form(client):
    with dataset_path('MOL-002').open('rb') as f:
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
        '/api/v1/chemspectra/predict/nmr_peaks_form',
        content_type='multipart/form-data',
        data=data
    )
    assert response.status_code == 200
    assert response.json['outline']['code'] == 200
    assert response.mimetype == 'application/json'


@mock.patch('requests.post', mock.Mock(return_value=ResponsePredictIr()))
def test_api_chemspectra_predict_infrared(client):
    with dataset_path('IR-004').open('rb') as f:
        spectrum_content = f.read()
    with dataset_path('MOL-002').open('rb') as f:
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
    assert response.json['outline']['code'] == 200
    assert response.mimetype == 'application/json'
