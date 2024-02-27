import json
from flask import (
    Blueprint, request, jsonify, abort,
)

# from chem_spectra.controller.helper.settings import get_ip_white_list
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.controller.helper.share import parse_array_to_dict_xys

from chem_spectra.model.inferencer import InferencerModel as InferModel
from chem_spectra.model.molecule import MoleculeModel
from chem_spectra.model.transformer import TransformerModel as TraModel

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
    mm = MoleculeModel(molfile, layout, decorate=True)

    outcome = InferModel.predict_nmr(
        mm=mm,
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
    layout = request.form.get('layout', default=None)
    peaks = request.form.get('peaks', default='{}')
    peaks = json.loads(peaks)
    shift = request.form.get('shift', default='{}')
    shift = json.loads(shift)
    molfile = FileContainer(request.files['molfile'])
    mm = MoleculeModel(molfile, layout, decorate=True)

    if not peaks:
        spectrum = FileContainer(request.files['spectrum'])
        if spectrum and layout == '13C':
            cv = TraModel(spectrum, molfile=molfile, params={'ext': 'jdx'}).to_converter()
            peaks = parse_array_to_dict_xys(cv.edit_peaks)

    if (not peaks) or (not molfile):
        abort(400)

    outcome = InferModel.predict_nmr(
        mm=mm,
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
    layout = request.form.get('layout', default=None)
    spectrum = FileContainer(request.files['spectrum'])
    molfile = FileContainer(request.files['molfile'])
    mm = MoleculeModel(molfile, layout)

    outcome = InferModel.predict_ir(
        mm=mm,
        spectrum=spectrum
    )
    if outcome:
        return jsonify(outcome)
    abort(400)


@infer_api.route('/predict/ms', methods=['POST'])
@infer_api.route('/api/v1/chemspectra/predict/ms', methods=['POST'])
def chemspectra_predict_ms():
    layout = request.form.get('layout', default=None)
    spectrum = FileContainer(request.files['spectrum'])
    molfile = FileContainer(request.files['molfile'])
    mm = MoleculeModel(molfile, layout)
    tm, _ = TraModel(spectrum, molfile=None, params={'ext': 'jdx'}).to_composer()
    if ((type(tm) is dict) and "invalid_molfile" in tm):
        return json.dumps(tm)

    outcome = InferModel.predict_ms(
        mm=mm,
        tm=tm
    )
    if outcome:
        return jsonify(outcome)
    abort(400)
