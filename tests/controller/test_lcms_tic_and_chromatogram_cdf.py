"""Two uploads that answered 500 on /file/convert.

- An LC/MS JCAMP carrying only a TIC converts, but the LC/MS preview draws
  the UV/Vis and MS traces and makes no image from a TIC alone.
  /file/convert then called `.read()` on `None`.
- An AIA chromatography netCDF (ANDI/CHROM, e.g. a diode-array trace) shares
  the .cdf extension with ANDI/MS but has no mass axis. Reading it raised
  KeyError: 'mass_values' on every endpoint.

The fixture and the netCDF are synthetic.
"""
import base64
import io

import netCDF4
import numpy as np
import pytest

TIC_ONLY = './tests/fixtures/source/lcms/lcms_tic_only.jdx'
ENDPOINTS = ('/api/v1/chemspectra/file/convert', '/zip_jcamp_n_img')


def _post(client, endpoint, content, name):
    return client.post(
        endpoint,
        content_type='multipart/form-data',
        data={'file': (io.BytesIO(content), name)},
    )


def test_a_tic_only_upload_converts_without_an_image(client):
    with open(TIC_ONLY, 'rb') as f:
        content = f.read()

    response = _post(client, ENDPOINTS[0], content, 'lcms_tic_only.jdx')

    assert response.status_code == 200
    body = response.get_json()
    assert body['img'] is None
    jcamp = base64.b64decode(body['jcamp']).decode('utf-8')
    assert '##DATA TYPE=LC/MS' in jcamp


def test_a_tic_only_upload_still_converts_on_zip_jcamp_n_img(client):
    with open(TIC_ONLY, 'rb') as f:
        content = f.read()

    response = _post(client, ENDPOINTS[1], content, 'lcms_tic_only.jdx')

    assert response.status_code == 200


def _chromatogram_cdf(path):
    with netCDF4.Dataset(path, 'w') as cdf:
        cdf.dataset_completeness = 'C1+C2'
        cdf.aia_template_revision = '1.0'
        cdf.detector_name = 'DAD1A,Sig=254,4  Ref=off'
        cdf.createDimension('point_number', 5)
        values = cdf.createVariable('ordinate_values', 'f4', ('point_number',))
        values[:] = np.array([0.0, 1.5, 4.0, 1.5, 0.0], dtype='f4')


@pytest.mark.parametrize('endpoint', ENDPOINTS)
def test_a_chromatogram_cdf_is_refused_with_the_reason(
        client, tmp_path, endpoint):
    path = tmp_path / 'dad.cdf'
    _chromatogram_cdf(path)

    response = _post(client, endpoint, path.read_bytes(), 'dad.cdf')

    assert response.status_code == 422
    error = response.get_json()['error']
    assert 'chromatogram' in error
    assert 'DAD1A' in error
    assert 'ANDI/MS' in error
