import io
import json


target_dir = './tests/fixtures/'
source_dir = 'source/'
file_jdx = '13C-DEPT135.dx'
bagit_file = 'cyclicvoltammetry/File053_BagIt.zip'
cv_file_jdx = 'cyclicvoltammetry/RCV_LSH-R444_full+Fc.jdx'

def test_api_chemspectra_file_convert_without_file(client):
    data = {}
    response = client.post(
        '/api/v1/chemspectra/file/convert',
        content_type='multipart/form-data',
        data=data
    )

    assert response.status_code == 400


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

def test_api_chemspectra_file_convert_bagit_format(client):
    with open(target_dir + source_dir + bagit_file, 'rb') as f:
        file_content = f.read()
    data = dict(
        file=(io.BytesIO(file_content), 'File053_BagIt.zip'),
    )
    response = client.post(
        '/api/v1/chemspectra/file/convert',
        content_type='multipart/form-data',
        data=data
    )

    assert response.status_code == 200
    assert response.mimetype == 'application/json'

    response_data_as_dict = json.loads(response.data)
    assert len(response_data_as_dict["list_jcamps"]) == 3
    

def test_api_chemspectra_file_save_single_file(client):
    with open(target_dir + source_dir + file_jdx, 'rb') as f:
        file_content = f.read()
    data = dict(
        src=(io.BytesIO(file_content), '13C-DEPT135.dx'),
        dst=(io.BytesIO(file_content), '13C-DEPT135.dx'),
        molfile=None
    )
    response = client.post(
        '/api/v1/chemspectra/file/save',
        content_type='multipart/form-data',
        data=data
    )

    assert response.status_code == 200
    assert response.mimetype == 'application/zip'

def test_api_chemspectra_file_save_multiple_files(client):
    with open(target_dir + source_dir + cv_file_jdx, 'rb') as f:
        file_content = f.read()
    data = dict(
        src=(io.BytesIO(file_content), 'cyclicvoltammetry/RCV_LSH-R444_full+Fc.jdx'),
        dst_list=[(io.BytesIO(file_content), 'cyclicvoltammetry/RCV_LSH-R444_full+Fc.jdx'), (io.BytesIO(file_content), 'cyclicvoltammetry/RCV_LSH-R444_full+Fc.jdx')],
        molfile=None
    )
    response = client.post(
        '/api/v1/chemspectra/file/save',
        content_type='multipart/form-data',
        data=data
    )

    assert response.status_code == 200
    assert response.mimetype == 'application/zip'

def test_api_chemspectra_file_refresh_single_file(client):
    with open(target_dir + source_dir + file_jdx, 'rb') as f:
        file_content = f.read()
    data = dict(
        dst=(io.BytesIO(file_content), '13C-DEPT135.dx'),
        molfile=None
    )
    response = client.post(
        '/api/v1/chemspectra/file/refresh',
        content_type='multipart/form-data',
        data=data
    )

    assert response.status_code == 200
    assert response.mimetype == 'application/json'

def test_api_chemspectra_file_refresh_multiple_files(client):
    with open(target_dir + source_dir + cv_file_jdx, 'rb') as f:
        file_content = f.read()
    data = dict(
        dst_list=[(io.BytesIO(file_content), 'cyclicvoltammetry/RCV_LSH-R444_full+Fc.jdx'), (io.BytesIO(file_content), 'cyclicvoltammetry/RCV_LSH-R444_full+Fc.jdx')],
        molfile=None
    )
    response = client.post(
        '/api/v1/chemspectra/file/refresh',
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
