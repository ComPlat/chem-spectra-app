import json
from flask import (
    Flask, Blueprint, request, redirect, jsonify, send_file, abort,
)

from chem_spectra.controller.helper.settings import get_ip_white_list
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.model.inferencer import InferencerModel as InfModel


infer_api = Blueprint('inference_api', __name__)


@infer_api.route('/predict/by_peaks_json', methods=['POST'])
@infer_api.route('/api/v1/chemspectra/predict/by_peaks_json', methods=['POST'])
def chemspectra_predict_by_peaks_json():
    try:
        payload = request.json
        layout = payload.get('layout')
        peaks = payload.get('peaks')
        shift = payload.get('shift')
        molfile = FileContainer().from_str(payload.get('molfile'))

        rsp = InfModel(layout, molfile, peaks, shift).predict_by_peaks()
        if rsp:
            return jsonify(
                status=True,
                result=rsp.json(),
            )
            abort(400)
    except:
        abort(500)


@infer_api.route('/predict/by_peaks_form', methods=['POST'])
@infer_api.route('/api/v1/chemspectra/predict/by_peaks_form', methods=['POST'])
def chemspectra_predict_by_peaks_form():
    try:
        molfile = FileContainer(request.files['molfile'])
        layout = request.form.get('layout', default=None)
        peaks = request.form.get('peaks', default=None)
        peaks = json.loads(peaks)
        shift = request.form.get('shift', default=None)
        shift = json.loads(shift)

        if (not peaks) or (not molfile):
            abort(400)

        rsp = InfModel(layout, molfile, peaks, shift).predict_by_peaks()
        if rsp:
            return jsonify(
                status=True,
                result=rsp.json(),
            )
            abort(400)
    except:
        abort(500)



@infer_api.route('/predict/infrared', methods=['POST'])
@infer_api.route('/api/v1/chemspectra/predict/infrared', methods=['POST'])
def chemspectra_predict_infrared():
    try:
        molfile = FileContainer(request.files['molfile'])
        spectrum = FileContainer(request.files['spectrum'])

        # results = InfraredModel(molfile, spectrum).check_fgs()

        # return jsonify(
        #     status=True,
        #     results=results,
        # )
        abort(400)
    except:
        abort(500)
