import json
from flask import jsonify, Blueprint, request
import os
import shutil
from flask_jwt_extended import jwt_required

spectra_layout_api = Blueprint('spectra_layout_api', __name__)
script_dir = os.path.dirname(__file__)
data_type_json_path = os.path.join(script_dir, '../lib/converter/jcamp/data_type.json')

def load_data_types():
    try:
        with open(data_type_json_path, 'r') as mapping_file:
            return json.load(mapping_file)
    except FileNotFoundError:
        example_json_path = os.path.join(script_dir, '../lib/converter/jcamp/data_type.json.example')
        shutil.copy(example_json_path, data_type_json_path)
        with open(data_type_json_path, 'r') as mapping_file:
            return json.load(mapping_file)

def save_data_types(data_types):
    with open(data_type_json_path, 'w') as mapping_file:
        json.dump(data_types, mapping_file, indent=4)

@spectra_layout_api.route('/api/v1/chemspectra/spectra_layouts', methods=['GET', 'PUT', 'DELETE'])
@jwt_required()
def update_or_fetch_mapping():
    if request.method == 'GET':
         existing_data_types = load_data_types()
         return jsonify(existing_data_types["datatypes"]), 200   
    
    elif request.method == 'PUT':
        request_data = request.get_json()
        new_data_type_mapping = request_data.get("new_data_type")
        existing_data_types = load_data_types()

        for layout, data_type in new_data_type_mapping.items():
            if data_type == '':
                return jsonify({"message": "Invalid Data Type"}), 400
            elif layout in existing_data_types["datatypes"] and data_type not in existing_data_types['datatypes'][layout]:
                existing_data_types["datatypes"][layout].append(data_type)
            elif layout in existing_data_types["datatypes"] and data_type in existing_data_types['datatypes'][layout]:
                return jsonify({"message": f"Data type '{data_type}' already exists"}), 400
            else:
                return jsonify({"message": f"Layout '{layout}' does not exist"}), 400
        
        save_data_types(existing_data_types)
        return jsonify({"message": "Data type created successfully"}), 200
    
    elif request.method == 'DELETE':
        request_data = request.get_json()
        data_type_mapping = request_data.get("data_type")
        existing_data_types = load_data_types()
        for layout, data_type in data_type_mapping.items():
            if layout in existing_data_types["datatypes"]:
                if data_type in existing_data_types['datatypes'][layout]:
                    existing_data_types['datatypes'][layout].remove(data_type)
                    save_data_types(existing_data_types)
                    return jsonify({"message": f"Data type '{data_type}' deleted successfully"}), 200
                else:
                    return jsonify({"message": f"Data type '{data_type}' not found in layout '{layout}'"}), 404
            else:
                return jsonify({"message": f"Layout '{layout}' does not exist"}), 400
    else:
        return jsonify({"message": "Method not allowed"}), 405
