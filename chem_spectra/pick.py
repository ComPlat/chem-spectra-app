from flask import (
    Flask, Blueprint, request, redirect, jsonify, send_file, abort
)
from .lib.spectra.helper import (
    allowed_file, convert2jcamp_img, convert2jcamp, convert2img
)
from .lib.process import to_zip_response
from .settings import get_ip_white_list


pk = Blueprint('pick', __name__)


@pk.before_app_request
def filter_remote_ip():
    trusted_servers = get_ip_white_list()
    if request.remote_addr not in trusted_servers:
        abort(403)


@pk.route('/peak_zip_jcamp_n_img', methods=['POST'])
def peak_zip_jcamp_n_img():
    try:
        file = request.files['file']
        if file: # and allowed_file(file):
            tf_jcamp, tf_img = convert2jcamp_img(file)
            memory = to_zip_response([tf_jcamp, tf_img])
            return send_file(
                memory,
                attachment_filename='spectrum.zip',
                as_attachment=True
            )
        abort(400)
    except:
        abort(500)


@pk.route('/edit_zip_jcamp_n_img', methods=['POST'])
def edit_zip_jcamp_n_img():
    try:
        file = request.files['file']
        peaks_str = request.form['peaks_str']
        if file: # and allowed_file(file):
            tf_jcamp, tf_img = convert2jcamp_img(file, peaks_str)
            memory = to_zip_response([tf_jcamp, tf_img])
            return send_file(
                memory,
                attachment_filename='spectrum.zip',
                as_attachment=True
            )
        abort(400)
    except:
        abort(500)


@pk.route('/peak_in_jcamp', methods=['POST'])
def peak_in_jcamp():
    try:
        file = request.files['file']
        if file: # and allowed_file(file):
            tf_jcamp = convert2jcamp(file)
            memory = to_zip_response([tf_jcamp])
            return send_file(
                memory,
                attachment_filename='spectrum.zip',
                as_attachment=True
            )
        abort(400)
    except:
        abort(500)


@pk.route('/peak_in_image', methods=['POST'])
def peak_in_image():
    try:
        file = request.files['file']
        if file: # and allowed_file(file):
            tf_img = convert2img(file)
            memory = to_zip_response([tf_img])
            return send_file(
                memory,
                attachment_filename='spectrum.zip',
                as_attachment=True
            )
        abort(400)
    except:
        abort(500)