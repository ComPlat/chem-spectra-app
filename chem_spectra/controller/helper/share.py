import io
import zipfile
import numpy as np
import json
from os.path import basename


ALLOWED_EXTENSIONS = set(['dx', 'jdx', 'raw', 'mzml'])


def allowed_file(file):
    ext = file.filename.split('.')[-1].lower()
    return ext in ALLOWED_EXTENSIONS


def get_fname(abs_path, filename):
    if not filename:
        return basename(abs_path)

    ext = abs_path.split('.')[-1]

    return filename + '.' + ext


def to_zip_response(tmp_arr, filename=False):
    memory = io.BytesIO()
    with zipfile.ZipFile(memory, 'w') as zf:
        for tmp in tmp_arr:
            abs_path = tmp.name
            fname = get_fname(abs_path, filename).replace(' ', '_')
            zf.write(abs_path, fname)
    memory.seek(0)
    for tmp in tmp_arr:
        tmp.close()
    return memory


def extract_params(request):
    scan = float(request.form.get('scan', default=0))
    scan = 0 if np.isnan(scan) else int(scan)
    thres = float(request.form.get('thres', default=0))
    mass = float(request.form.get('mass', default=0))
    clear = bool(request.form.get('clear', default=False))
    predict = request.form.get('predict', default='{}')

    params = {
        'peaks_str': request.form.get('peaks_str', default=None),
        'select_x': request.form.get('shift_select_x', default=None),
        'ref_name': request.form.get('shift_ref_name', default=None),
        'ref_value': request.form.get('shift_ref_value', default=None),
        'scan': scan,
        'thres': thres,
        'mass': mass,
        'molfile': request.form.get('molfile', default=None),
        'clear': clear,
        'predict': predict,
    }
    has_params = (
        params.get('peaks_str') or
        params.get('select_x') or
        params.get('ref_name') or
        params.get('ref_value') or
        params.get('scan') or
        params.get('thres') or
        params.get('mass') or
        params.get('molfile') or
        params.get('clear') or
        params.get('predict')
    )
    if not has_params:
        params = False
    return params
