from flask import (
    Flask, Blueprint, request, redirect, jsonify, send_file, abort, jsonify
)
import base64
from .lib.spectra.helper import (
    allowed_file, convert2jcamp_img, convert2jcamp, convert2img
)
from .lib.process import to_zip_response, extract_params
from .settings import get_ip_white_list


pk = Blueprint('pick', __name__)


@pk.before_app_request
def filter_remote_ip():
    trusted_servers = get_ip_white_list()
    if request.remote_addr not in trusted_servers:
        abort(403)


@pk.route('/zip_jcamp_n_img', methods=['POST'])
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


@pk.route('/zip_jcamp', methods=['POST'])
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


@pk.route('/zip_image', methods=['POST'])
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


@pk.route('/jcamp', methods=['POST'])
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


@pk.route('/image', methods=['POST'])
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


@pk.route('/api/v1/chemspectra/file/convert', methods=['POST'])
def chemspectra_file_convert():
    try:
        file = request.files['file']
        params = extract_params(request)
        if file: # and allowed_file(file):
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
        return jsonify(
            status=False,
            jcamp='',
            img=''
        )
        abort(500)


@pk.route('/api/v1/chemspectra/file/save', methods=['POST'])
def chemspectra_file_save():
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
        return jsonify(
            status=False,
            jcamp='',
            img=''
        )
        abort(500)
