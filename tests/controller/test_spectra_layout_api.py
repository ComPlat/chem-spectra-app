import json
from unittest.mock import patch

test_json_path = './tests/fixtures/test_data_types.json'
orig_json_path = 'chem_spectra.controller.spectra_layout_api.data_type_json_path'

def test_fetch_mapping_get_spectra_layouts(client):

    with patch(orig_json_path, new=test_json_path):
        response = client.get('/api/v1/chemspectra/spectra_layouts')
        response_data = response.json
        assert response.status_code == 200
        assert response_data['datatypes'] == {
            "INFRARED": ["INFRARED SPECTRUM"],
            "MS": ["MASS SPECTRUM"],
            "NMR": ["NMR SPECTRUM", "NMRSPECTRUM"],
            "RAMAN": ["RAMAN"]
        }

def test_create_data_type(client):
    access_token = fetch_access_token(client)
    new_data_type = {
        "new_data_type": {
            "MS": "MASS SPEC"
        }
    }

    with patch(orig_json_path, test_json_path):
        response = client.put('/api/v1/chemspectra/spectra_layouts', json=new_data_type, 
                               headers={'Authorization': f'Bearer {access_token}'})

        assert response.status_code == 200
        
        response_data = response.json
        assert "message" in response_data
        response_data["message"] == "Data type created successfully"
        
        delete_data_type = {
        "data_type": {
            "MS": "MASS SPEC"
        }
    }
        response = client.delete('/api/v1/chemspectra/spectra_layouts', json=delete_data_type,
                                 headers={'Authorization': f'Bearer {access_token}'})

def test_create_data_type_unchanged(client): 
    access_token = fetch_access_token(client)   
    # data type already exists in JSON
    new_data_type = {
        "new_data_type": {
             "INFRARED": "INFRARED SPECTRUM"
        }
    }

    with patch(orig_json_path, test_json_path):
        response = client.put('/api/v1/chemspectra/spectra_layouts', json=new_data_type, 
                               headers={'Authorization': f'Bearer {access_token}'})
        assert response.status_code == 400

        response_data = response.json
        
        assert "message" in response_data
        assert response_data["message"] == "Data type 'INFRARED SPECTRUM' already exists"

def test_create_data_type_layout_does_not_exist(client):
    access_token = fetch_access_token(client)
    new_data_type = {
        "new_data_type": {
             "UNKNOWN": "NA"
        }
    }

    with patch(orig_json_path, test_json_path):
        response = client.put('/api/v1/chemspectra/spectra_layouts', json=new_data_type, 
                               headers={'Authorization': f'Bearer {access_token}'})
        assert response.status_code == 400

        response_data = response.json
        
        assert "message" in response_data
        assert response_data["message"] == "Layout 'UNKNOWN' does not exist"

def test_delete_data_type(client):
    access_token = fetch_access_token(client)
    # create a new data type
    new_data_type = {
        "new_data_type": {
            "INFRARED": "IR"
        }
    }
    with patch(orig_json_path, test_json_path):
        response = client.put('/api/v1/chemspectra/spectra_layouts', json=new_data_type, 
                               headers={'Authorization': f'Bearer {access_token}'})
        assert response.status_code == 200

        # delete the new data type
        delete_data_type = {
            "data_type": {
            "INFRARED": "IR"
            }
        }   

        response = client.delete('/api/v1/chemspectra/spectra_layouts', json=delete_data_type,
                                 headers={'Authorization': f'Bearer {access_token}'})
        assert response.status_code == 200
        response_data = response.json

        assert "message" in response_data
        assert response_data["message"] == "Data type 'IR' deleted successfully"
