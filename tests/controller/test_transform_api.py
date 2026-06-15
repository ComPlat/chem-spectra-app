import io
import json
import zipfile

from tests.dataset_catalog import dataset_path, legacy_source_path

file_jdx = '13C-DEPT135.dx'
result_dir = 'result/'

peaks_str = '745.0957757310398,0.2787140606224312#1018.4864309069585,0.31625977127489585#1154.473492548866,0.32047998816450246'  # noqa:


def test_zip_jcamp_n_img(client):
    with open(legacy_source_path(file_jdx), 'rb') as f:
        file_content = f.read()
    data = dict(
        file=(io.BytesIO(file_content), '13C-DEPT135.dx'),
    )
    response = client.post(
        '/zip_jcamp_n_img',
        content_type='multipart/form-data',
        data=data
    )

    assert response.status_code == 200
    assert response.mimetype == 'application/zip'

    header_json = json.loads(response.headers['X-Extra-Info-JSON'])
    assert header_json['invalid_molfile'] is False

def test_zip_jcamp_n_img_invalid_molfile(client):
    with open(legacy_source_path(file_jdx), 'rb') as f:
        file_content = f.read()

    with dataset_path('MOL-001').open('rb') as f:
        molfile = f.read()

    data = dict(
        file=(io.BytesIO(file_content), '13C-DEPT135.dx'),
        molfile=(io.BytesIO(molfile), 'invalid_molfile.mol'),
        simulatenmr='true'
    )
    
    response = client.post(
        '/zip_jcamp_n_img',
        content_type='multipart/form-data',
        data=data
    )

    assert response.status_code == 200
    assert response.mimetype == 'application/zip'

    header_json = json.loads(response.headers['X-Extra-Info-JSON'])
    assert header_json['invalid_molfile'] is True


def test_zip_jcamp_n_img_with_peaks_str(client):
    with open(legacy_source_path(file_jdx), 'rb') as f:
        file_content = f.read()
    data = dict(
        file=(io.BytesIO(file_content), '13C-DEPT135.dx'),
        peaks_str=peaks_str
    )
    response = client.post(
        '/zip_jcamp_n_img',
        content_type='multipart/form-data',
        data=data
    )

    assert response.status_code == 200
    assert response.mimetype == 'application/zip'


def test_zip_jcamp(client):
    with open(legacy_source_path(file_jdx), 'rb') as f:
        file_content = f.read()
    data = dict(
        file=(io.BytesIO(file_content), '13C-DEPT135.dx'),
    )
    response = client.post(
        '/zip_jcamp',
        content_type='multipart/form-data',
        data=data
    )

    assert response.status_code == 200
    assert response.mimetype == 'application/zip'


def test_zip_image(client):
    with open(legacy_source_path(file_jdx), 'rb') as f:
        file_content = f.read()
    data = dict(
        file=(io.BytesIO(file_content), '13C-DEPT135.dx'),
    )
    response = client.post(
        '/zip_image',
        content_type='multipart/form-data',
        data=data
    )

    assert response.status_code == 200
    assert response.mimetype == 'application/zip'
    
def test_combine_images_no_file(client):
    response = client.post(
        '/combine_images',
        content_type='multipart/form-data',
        data=None
    )
    
    assert response.status_code == 400

def test_combine_images(client):
    #TODO: implement later
    pass


def test_image(client):
    with open(legacy_source_path(file_jdx), 'rb') as f:
        file_content = f.read()
    data = dict(
        file=(io.BytesIO(file_content), '13C-DEPT135.dx'),
    )
    response = client.post(
        '/image',
        content_type='multipart/form-data',
        data=data
    )

    assert response.status_code == 200
    assert response.mimetype == 'image/png'
