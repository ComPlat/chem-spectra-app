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
