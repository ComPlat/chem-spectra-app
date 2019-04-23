import pytest
import io


target_dir = './tests/fixtures/'
source_dir = 'source/'
file_jdx = '13C-DEPT135.dx'

peaks_str = '745.0957757310398,0.2787140606224312#1018.4864309069585,0.31625977127489585#1154.473492548866,0.32047998816450246'


def test_zip_jcamp_n_img(client):
    with open(target_dir + source_dir + file_jdx, 'rb') as f:
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


def test_zip_jcamp_n_img_with_peaks_str(client):
    with open(target_dir + source_dir + file_jdx, 'rb') as f:
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
    with open(target_dir + source_dir + file_jdx, 'rb') as f:
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
    with open(target_dir + source_dir + file_jdx, 'rb') as f:
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


def test_jcamp(client):
    with open(target_dir + source_dir + file_jdx, 'rb') as f:
        file_content = f.read()
    data = dict(
        file=(io.BytesIO(file_content), '13C-DEPT135.dx'),
    )
    response = client.post(
        '/jcamp',
        content_type='multipart/form-data',
        data=data
    )

    assert response.status_code == 200
    assert response.mimetype == 'application/octet-stream'


def test_image(client):
    with open(target_dir + source_dir + file_jdx, 'rb') as f:
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
