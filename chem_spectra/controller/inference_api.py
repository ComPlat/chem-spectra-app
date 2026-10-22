import json
from flask import (
    Blueprint, request, jsonify,
)

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.controller.helper.refusal import uploaded_file
from chem_spectra.controller.helper.share import parse_array_to_dict_xys

from chem_spectra.lib.converter.ms import MSConversionFailed
from chem_spectra.lib.converter.share import UnconvertibleSpectrum
from chem_spectra.model.inferencer import (
    InferencerModel as InferModel, outline, UNUSABLE_REQUEST, NOT_AVAILABLE,
)
from chem_spectra.model.molecule import MoleculeModel, UnreadableMolfile
from chem_spectra.model.transformer import TransformerModel as TraModel

infer_api = Blueprint('inference_api', __name__)

# Every answer from these routes, result or failure, carries `outline`, and
# goes out through `_respond`.


def _respond(result):
    # HTTP 200 whatever `outline.code` says. The ELN discards any other status
    # (`rsp.code == 200 ? rsp.parsed_response : nil`), so a failure sent as,
    # say, 503 would reach the user as nothing at all. Once the ELN reads the
    # body of a non-200 answer, return `result['outline']['code']` here
    # instead; nothing else needs to change.
    return jsonify(result), 200


def _unusable(text):
    return _respond(outline(UNUSABLE_REQUEST, text))


def _molecule(molfile, layout, decorate=False):
    """The molecule, or None when there is no readable molfile."""
    if molfile is None or not molfile.core:
        return None
    try:
        return MoleculeModel(molfile, layout, decorate=decorate)
    except UnreadableMolfile:
        return None


def _uploaded(field):
    uploaded = uploaded_file(request, field)
    return FileContainer(uploaded) if uploaded else None


_UNPARSABLE = object()


def _json_field(name, default):
    """The form field parsed as JSON; `default` when absent or `null`."""
    try:
        value = json.loads(request.form.get(name, default=default))
    except ValueError:
        return _UNPARSABLE
    return json.loads(default) if value is None else value


def _peaks_problem(peaks):
    """Why `peaks` cannot be sent to a predictor, or None."""
    if not peaks:
        return 'No peaks were sent.'
    if not isinstance(peaks, list):
        return 'peaks must be a list of {"x": ..., "y": ...} objects.'
    for peak in peaks:
        x = peak.get('x') if isinstance(peak, dict) else None
        if isinstance(x, bool) or not isinstance(x, (int, float)):
            return 'Every peak needs a numeric "x".'
    return None


def _conversion_failed(err):
    """The outline for a spectrum that could not be converted."""
    if isinstance(err, MSConversionFailed) and \
            getattr(err, 'status', 422) >= 500:
        return _respond(outline(NOT_AVAILABLE, str(err)))
    return _unusable('The spectrum could not be read: {}'.format(err))


def _predict_nmr(mm, layout, peaks, shift):
    problem = _peaks_problem(peaks)
    if problem:
        return _unusable(problem)
    if not isinstance(shift, dict):
        return _unusable('shift must be a JSON object.')
    return _respond(InferModel.predict_nmr(
        mm=mm,
        layout=layout,
        peaks=peaks,
        shift=shift
    ))


@infer_api.route('/predict/by_peaks_json', methods=['POST'])
@infer_api.route(
    '/api/v1/chemspectra/predict/nmr_peaks_json', methods=['POST']
)
def chemspectra_predict_by_peaks_json():
    payload = request.get_json(silent=True)
    if not isinstance(payload, dict):
        return _unusable('The request body is not a JSON object.')
    layout = payload.get('layout')
    peaks = payload.get('peaks')
    shift = payload.get('shift')
    shift = {} if shift is None else shift
    moltxt = payload.get('molfile')
    molfile = FileContainer().from_str(moltxt) \
        if isinstance(moltxt, str) and moltxt else None
    mm = _molecule(molfile, layout, decorate=True)
    if not mm:
        return _unusable('No readable molfile was sent.')

    return _predict_nmr(mm, layout, peaks, shift)


@infer_api.route('/predict/by_peaks_form', methods=['POST'])
@infer_api.route(
    '/api/v1/chemspectra/predict/nmr_peaks_form', methods=['POST']
)
def chemspectra_predict_by_peaks_form():
    layout = request.form.get('layout', default=None)
    peaks = _json_field('peaks', '{}')
    shift = _json_field('shift', '{}')
    if peaks is _UNPARSABLE or shift is _UNPARSABLE:
        return _unusable('peaks and shift must be JSON.')
    molfile = _uploaded('molfile')
    mm = _molecule(molfile, layout, decorate=True)
    if not mm:
        return _unusable('No readable molfile was sent.')

    if not peaks:
        spectrum = _uploaded('spectrum')
        if spectrum and layout == '13C':
            try:
                cv = TraModel(
                    spectrum, molfile=molfile, params={'ext': 'jdx'}
                ).to_converter()
            except (UnconvertibleSpectrum, MSConversionFailed) as err:
                return _conversion_failed(err)
            if not cv:
                return _unusable('The spectrum could not be read.')
            peaks = parse_array_to_dict_xys(cv.edit_peaks)

    if not peaks:
        return _unusable(
            'No peaks were sent, and none could be picked from a 13C '
            'spectrum.'
        )

    return _predict_nmr(mm, layout, peaks, shift)


@infer_api.route('/predict/infrared', methods=['POST'])
@infer_api.route('/api/v1/chemspectra/predict/infrared', methods=['POST'])
def chemspectra_predict_infrared():
    layout = request.form.get('layout', default=None)
    spectrum = _uploaded('spectrum')
    if not spectrum:
        return _unusable('No spectrum was sent.')
    mm = _molecule(_uploaded('molfile'), layout)
    if not mm:
        return _unusable('No readable molfile was sent.')

    return _respond(InferModel.predict_ir(
        mm=mm,
        spectrum=spectrum
    ))


@infer_api.route('/predict/ms', methods=['POST'])
@infer_api.route('/api/v1/chemspectra/predict/ms', methods=['POST'])
def chemspectra_predict_ms():
    layout = request.form.get('layout', default=None)
    spectrum = _uploaded('spectrum')
    if not spectrum:
        return _unusable('No spectrum was sent.')
    mm = _molecule(_uploaded('molfile'), layout)
    if not mm:
        return _unusable('No readable molfile was sent.')
    try:
        tm, _ = TraModel(
            spectrum, molfile=None, params={'ext': 'jdx'}
        ).to_composer()
    except (UnconvertibleSpectrum, MSConversionFailed) as err:
        return _conversion_failed(err)
    if not tm:
        return _unusable('The spectrum could not be read.')
    if not hasattr(tm, 'prism_peaks'):
        # a JCAMP that converts, but not to a mass spectrum
        return _unusable('MS prediction needs a mass spectrum.')

    return _respond(InferModel.predict_ms(
        mm=mm,
        tm=tm
    ))
