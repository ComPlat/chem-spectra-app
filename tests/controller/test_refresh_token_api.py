def test_login(client):
    response = client.post('/api/v1/chemspectra/login')
    data = response.json
    assert response.status_code == 200
    assert "access_token" in data
    assert "refresh_token" in data

def test_refresh(client):
    # get refresh token
    response = client.post('/api/v1/chemspectra/login')
    data = response.json
    refresh_token = data.get('refresh_token')

    refresh_response = client.post('/api/v1/chemspectra/refresh', headers={'Authorization': f'Bearer {refresh_token}'})
    data = refresh_response.json
    assert "access_token" in data
    assert refresh_response.status_code == 200

def test_refresh_without_refresh_token(client):
    response = client.post('/api/v1/chemspectra/refresh')
    assert response.status_code == 401