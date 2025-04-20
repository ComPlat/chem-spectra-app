import json
from flask import jsonify, Blueprint, request
import os
import shutil

spectra_layout_api = Blueprint('spectra_layout_api', __name__)
script_dir = os.path.dirname(__file__)
data_type_json_path = os.path.join(
    script_dir, '../lib/converter/jcamp/data_type.json')


def load_data_types():
    try:
        with open(data_type_json_path, 'r') as mapping_file:
            return json.load(mapping_file)
    except FileNotFoundError:
        example_json_path = os.path.join(
            script_dir, '../lib/converter/jcamp/data_type.json.example')
        shutil.copy(example_json_path, data_type_json_path)
        with open(data_type_json_path, 'r') as mapping_file:
            return json.load(mapping_file)


@spectra_layout_api.route('/api/v1/chemspectra/spectra_layouts', methods=['GET'])
def fetch_mapping():
    if request.method == 'GET':
        existing_data_types = load_data_types()
        return jsonify(existing_data_types), 200
    else:
        return jsonify({"message": "Method not allowed"}), 405
