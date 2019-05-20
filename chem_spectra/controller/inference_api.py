import json
from flask import (
    Blueprint, request, jsonify, abort,
)

# from chem_spectra.controller.helper.settings import get_ip_white_list
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.model.inferencer import InferencerModel as InferModel

infer_api = Blueprint('inference_api', __name__)


@infer_api.route('/predict/by_peaks_json', methods=['POST'])
@infer_api.route(
    '/api/v1/chemspectra/predict/nmr_peaks_json', methods=['POST']
)
def chemspectra_predict_by_peaks_json():
    payload = request.json
    layout = payload.get('layout')
    peaks = payload.get('peaks')
    shift = payload.get('shift')
    molfile = FileContainer().from_str(payload.get('molfile'))

    outcome = InferModel.predict_nmr(
        molfile=molfile,
        layout=layout,
        peaks=peaks,
        shift=shift
    )
    if outcome:
        return jsonify(outcome)
    abort(400)


@infer_api.route('/predict/by_peaks_form', methods=['POST'])
@infer_api.route(
    '/api/v1/chemspectra/predict/nmr_peaks_form', methods=['POST']
)
def chemspectra_predict_by_peaks_form():
    molfile = FileContainer(request.files['molfile'])
    layout = request.form.get('layout', default=None)
    peaks = request.form.get('peaks', default='{}')
    peaks = json.loads(peaks)
    shift = request.form.get('shift', default='{}')
    shift = json.loads(shift)

    if (not peaks) or (not molfile):
        abort(400)

    outcome = InferModel.predict_nmr(
        molfile=molfile,
        layout=layout,
        peaks=peaks,
        shift=shift
    )
    if outcome:
        return jsonify(outcome)
    abort(400)


@infer_api.route('/predict/infrared', methods=['POST'])
@infer_api.route('/api/v1/chemspectra/predict/infrared', methods=['POST'])
def chemspectra_predict_infrared():
    molfile = FileContainer(request.files['molfile'])
    spectrum = FileContainer(request.files['spectrum'])
    outcome = InferModel.predict_ir(
        molfile=molfile,
        spectrum=spectrum
    )
    if outcome:
        return jsonify(outcome)
    abort(400)
