import base64
from flask import (
    Blueprint, request, jsonify, send_file, abort,
)

# from chem_spectra.controller.helper.settings import get_ip_white_list
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.controller.helper.share import (
    to_zip_response, extract_params
)
from chem_spectra.model.transformer import TransformerModel as TraModel
from chem_spectra.model.molecule import MoleculeModel


file_api = Blueprint('file_api', __name__)


@file_api.route('/api/v1/chemspectra/file/convert', methods=['POST'])
def chemspectra_file_convert():
    file = FileContainer(request.files['file'])
    params = extract_params(request)
    if file:
        tf_jcamp, tf_img = TraModel(file, params).convert2jcamp_img()
        if not tf_jcamp:
            abort(400)
        jcamp = base64.b64encode(tf_jcamp.read()).decode("utf-8")
        img = base64.b64encode(tf_img.read()).decode("utf-8")
        return jsonify(
            status=True,
            jcamp=jcamp,
            img=img
        )


@file_api.route('/api/v1/chemspectra/file/save', methods=['POST'])
def chemspectra_file_save():
    file = FileContainer(request.files['file'])
    filename = request.form.get('filename', default=None)
    params = extract_params(request)
    if file:  # and allowed_file(file):
        tm = TraModel(file, params)
        tf_jcamp, tf_img = tm.convert2jcamp_img()
        tf_arr = [tf_jcamp, tf_img, tm.tf_predict()]
        memory = to_zip_response(tf_arr, filename)
        return send_file(
            memory,
            attachment_filename='spectrum.zip',
            as_attachment=True
        )


@file_api.route('/api/v1/chemspectra/molfile/convert', methods=['POST'])
def chemspectra_molfile_convert():
    molfile = FileContainer(request.files['molfile'])
    mm = MoleculeModel(molfile)
    return jsonify(
        status=True,
        smi=mm.smi,
        mass=mm.mass
    )
