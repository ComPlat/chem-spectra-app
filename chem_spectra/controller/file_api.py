import json
import base64
from flask import (
    Flask, Blueprint, request, jsonify, send_file, abort,
)

from chem_spectra.controller.helper.settings import get_ip_white_list
from chem_spectra.controller.helper.share import (
    allowed_file, to_zip_response, extract_params
)

from chem_spectra.model.transformer import TransformerModel as TraModel
from chem_spectra.model.molecule import MoleculeModel


file_api = Blueprint('file_api', __name__)


@file_api.route('/api/v1/chemspectra/file/convert', methods=['POST'])
def chemspectra_file_convert():
    try:
        file = request.files['file']
        params = extract_params(request)
        if file:
            tf_jcamp, tf_img = TraModel(file, params).convert2jcamp_img()
            jcamp = base64.b64encode(tf_jcamp.read()).decode("utf-8")
            img = base64.b64encode(tf_img.read()).decode("utf-8")
            return jsonify(
                status=True,
                jcamp=jcamp,
                img=img
            )
        abort(400)
    except:
        abort(500)


@file_api.route('/api/v1/chemspectra/file/save', methods=['POST'])
def chemspectra_file_save():
    try:
        file = request.files['file']
        filename = request.form.get('filename', default=None)
        params = extract_params(request)
        if file: # and allowed_file(file):
            tf_jcamp, tf_img = TraModel(file, params).convert2jcamp_img()
            memory = to_zip_response([tf_jcamp, tf_img], filename)
            return send_file(
                memory,
                attachment_filename='spectrum.zip',
                as_attachment=True
            )
        abort(400)
    except:
        abort(500)


@file_api.route('/api/v1/chemspectra/molfile/convert', methods=['POST'])
def chemspectra_molfile_convert():
    try:
        molfile = request.files['molfile'].stream.read().decode('utf-8')
        mm = MoleculeModel(molfile)
        return jsonify(
            status=True,
            smi=mm.can(),
            mass=mm.mass()
        )
    except:
        abort(500)
