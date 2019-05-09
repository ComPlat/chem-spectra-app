import io
import zipfile
import json
import math
from os.path import basename


ALLOWED_EXTENSIONS = set(['dx', 'jdx', 'raw', 'mzml', 'jcamp'])


def allowed_file(file):
    ext = file.filename.split('.')[-1].lower()
    return ext in ALLOWED_EXTENSIONS


def get_fname(abs_path, filename):
    if not filename:
        return basename(abs_path)

    ext = abs_path.split('.')[-1]

    return filename + '.' + ext


def to_zip_response(src_tmp_arr, filename=False):
    tmp_arr = [el for el in src_tmp_arr if el]
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


def parse_float(val, default):
    try:
        return float(val)
    except:
        return float(default)


def extract_params(request):
    scan = parse_float(request.form.get('scan', default=0), 0)
    scan = 0 if math.isnan(scan) else int(scan)
    thres = parse_float(request.form.get('thres', default=0.0), 0.0)
    mass = parse_float(request.form.get('mass', default=1.0), 1.0)
    clear = bool(request.form.get('clear', default=False))
    ext = request.form.get('ext', default='')
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
        'ext': ext,
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
        params.get('predict') or
        params.get('ext')
    )
    if not has_params:
        params = False
    return params
