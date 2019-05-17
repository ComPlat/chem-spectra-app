import io


target_dir = './tests/fixtures/'
source_dir = 'source/'
file_jdx = '13C-DEPT135.dx'


def test_api_chemspectra_file_convert(client):
    with open(target_dir + source_dir + file_jdx, 'rb') as f:
        file_content = f.read()
    data = dict(
        file=(io.BytesIO(file_content), '13C-DEPT135.dx'),
    )
    response = client.post(
        '/api/v1/chemspectra/file/convert',
        content_type='multipart/form-data',
        data=data
    )

    assert response.status_code == 200
    assert response.mimetype == 'application/json'


def test_api_chemspectra_file_save(client):
    with open(target_dir + source_dir + file_jdx, 'rb') as f:
        file_content = f.read()
    data = dict(
        file=(io.BytesIO(file_content), '13C-DEPT135.dx'),
    )
    response = client.post(
        '/api/v1/chemspectra/file/save',
        content_type='multipart/form-data',
        data=data
    )

    assert response.status_code == 200
    assert response.mimetype == 'application/zip'


def test_api_chemspectra_molfile_convert(client):
    with open(target_dir + source_dir + '/molfile/svs813f1_B.mol', 'rb') as f:
        file_content = f.read()
    data = dict(
        molfile=(io.BytesIO(file_content), 'svs813f1_B.mol'),
    )
    response = client.post(
        '/api/v1/chemspectra/molfile/convert',
        content_type='multipart/form-data',
        data=data
    )
    assert response.status_code == 200
    assert response.mimetype == 'application/json'
