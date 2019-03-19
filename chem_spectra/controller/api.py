from flask import (
    Flask, Blueprint, request, redirect, jsonify, send_file, abort, jsonify
)
import base64
from chem_spectra.controller.helper import (
    allowed_file, convert2jcamp_img, convert2jcamp, convert2img,
)
from chem_spectra.controller.settings import get_ip_white_list
from chem_spectra.model.process import to_zip_response, extract_params
from chem_spectra.model.predict import predict_by_peaks
from chem_spectra.model.converter.chem import molfile2chem
import json

ctrl = Blueprint('api', __name__)


@ctrl.before_app_request
def filter_remote_ip():
    trusted_servers = get_ip_white_list()
    # if request.remote_addr not in trusted_servers:
    #     abort(403)


@ctrl.route('/zip_jcamp_n_img', methods=['POST'])
def zip_jcamp_n_img():
    try:
        file = request.files['file']
        params = extract_params(request)
        if file: # and allowed_file(file):
            tf_jcamp, tf_img = convert2jcamp_img(file, params)
            memory = to_zip_response([tf_jcamp, tf_img])
            return send_file(
                memory,
                attachment_filename='spectrum.zip',
                as_attachment=True
            )
        abort(400)
    except:
        abort(500)


@ctrl.route('/zip_jcamp', methods=['POST'])
def zip_jcamp():
    try:
        file = request.files['file']
        params = extract_params(request)
        if file: # and allowed_file(file):
            tf_jcamp = convert2jcamp(file, params)
            memory = to_zip_response([tf_jcamp])
            return send_file(
                memory,
                attachment_filename='spectrum.zip',
                as_attachment=True
            )
        abort(400)
    except:
        abort(500)


@ctrl.route('/zip_image', methods=['POST'])
def zip_image():
    try:
        file = request.files['file']
        params = extract_params(request)
        if file: # and allowed_file(file):
            tf_img = convert2img(file, params)
            memory = to_zip_response([tf_img])
            return send_file(
                memory,
                attachment_filename='spectrum.zip',
                as_attachment=True
            )
        abort(400)
    except:
        abort(500)


@ctrl.route('/jcamp', methods=['POST'])
def jcamp():
    try:
        file = request.files['file']
        params = extract_params(request)
        if file: # and allowed_file(file):
            tf_jcamp = convert2jcamp(file, params)
            return send_file(
                tf_jcamp,
                attachment_filename='spectrum.jdx',
                as_attachment=True
            )
        abort(400)
    except:
        abort(500)


@ctrl.route('/image', methods=['POST'])
def image():
    try:
        file = request.files['file']
        params = extract_params(request)
        if file: # and allowed_file(file):
            tf_img = convert2img(file, params)
            return send_file(
                tf_img,
                attachment_filename='spectrum.png',
                as_attachment=True,
                mimetype='image/png'
            )
        abort(400)
    except:
        abort(500)


@ctrl.route('/api/v1/chemspectra/file/convert', methods=['POST'])
def chemspectra_file_convert():
    try:
        file = request.files['file']
        params = extract_params(request)
        if file:
            tf_jcamp, tf_img = convert2jcamp_img(file, params)
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


@ctrl.route('/api/v1/chemspectra/file/save', methods=['POST'])
def chemspectra_file_save():
    try:
        file = request.files['file']
        filename = request.form.get('filename', default=None)
        params = extract_params(request)
        if file: # and allowed_file(file):
            tf_jcamp, tf_img = convert2jcamp_img(file, params)
            memory = to_zip_response([tf_jcamp, tf_img], filename)
            return send_file(
                memory,
                attachment_filename='spectrum.zip',
                as_attachment=True
            )
        abort(400)
    except:
        abort(500)


@ctrl.route('/predict/by_peaks_json', methods=['POST'])
@ctrl.route('/api/v1/chemspectra/predict/by_peaks_json', methods=['POST'])
def chemspectra_predict_by_peaks_json():
    try:
        payload = request.json
        layout = payload.get('layout')
        molfile = payload.get('molfile')
        peaks = payload.get('peaks')
        shift = payload.get('shift')

        rsp = predict_by_peaks(layout, molfile, peaks, shift)
        if rsp:
            return jsonify(
                status=True,
                result=rsp.json(),
            )
        abort(400)
    except:
        abort(500)


@ctrl.route('/predict/by_peaks_form', methods=['POST'])
@ctrl.route('/api/v1/chemspectra/predict/by_peaks_form', methods=['POST'])
def chemspectra_predict_by_peaks_form():
    try:
        molfile = request.files['molfile']
        molfile = molfile.stream.read().decode('utf-8')
        layout = request.form.get('layout', default=None)
        peaks = request.form.get('peaks', default=None)
        peaks = json.loads(peaks)
        shift = request.form.get('shift', default=None)
        shift = json.loads(shift)

        if (not peaks) or (not molfile):
            abort(400)

        rsp = predict_by_peaks(layout, molfile, peaks, shift)
        if rsp:
            return jsonify(
                status=True,
                result=rsp.json(),
            )
        abort(400)
    except:
        abort(500)


@ctrl.route('/api/v1/chemspectra/molfile/convert', methods=['POST'])
def chemspectra_molfile_convert():
    try:
        molfile = request.files['molfile']
        smi, mass = molfile2chem(molfile)
        return jsonify(
            status=True,
            smi=smi,
            mass=mass
        )
    except:
        abort(500)
