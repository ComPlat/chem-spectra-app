import io
import zipfile
import math
import os.path as os_path
from os.path import basename


ALLOWED_EXTENSIONS = set(['dx', 'jdx', 'raw', 'mzml', 'mzxml', 'jcamp'])


def allowed_file(file):
    ext = file.filename.split('.')[-1].lower()
    return ext in ALLOWED_EXTENSIONS


def get_fname(abs_path, filename, is_src):
    if not filename:
        return basename(abs_path)

    ext = abs_path.split('.')[-1]
    header = 'orig_' if is_src else ''

    return '{}{}.{}'.format(header, filename, ext)


def to_zip_response(src_tmp_arr, filename=False, src_idx=-1):
    tmp_arr = [el for el in src_tmp_arr if el]
    memory = io.BytesIO()
    with zipfile.ZipFile(memory, 'w') as zf:
        for idx, tmp in enumerate(tmp_arr):
            abs_path = tmp.name
            is_src = idx == src_idx
            fname = get_fname(abs_path, filename, is_src).replace(' ', '_')
            zf.write(abs_path, fname)
    memory.seek(0)
    for tmp in tmp_arr:
        tmp.close()
    return memory


def to_zip_bag_it_response(src_tmp_arr, filename=False, src_idx=-1):
    tmp_arr = [el for el in src_tmp_arr if el]
    memory = io.BytesIO()
    with zipfile.ZipFile(memory, 'w') as zf:
        for idx_sub, sub_arr in enumerate(tmp_arr):
            folder_path = 'curve_' + str(idx_sub)
            tmp_sub_arr = [el for el in sub_arr if el]
            for idx, tmp in enumerate(tmp_sub_arr):
                abs_path = tmp.name
                is_src = idx == src_idx
                fname = get_fname(abs_path, filename, is_src).replace(' ', '_')
                fname = os_path.join(folder_path, fname)
                zf.write(abs_path, fname)
            for tmp in tmp_sub_arr:
                tmp.close()
    memory.seek(0)

    return memory


def parse_float(val, default):
    try:
        return float(val)
    except:  # noqa
        return float(default)

def parse_int(val, default):
    try:
        return int(val)
    except:  # noqa
        return int(default)


def parse_fname(request):
    fil_name = request.files.get('file') and request.files.get('file').filename
    src_name = request.files.get('src') and request.files.get('src').filename
    def_name = request.form.get('fname', default=False)
    return def_name or fil_name or src_name or ''


def extract_params(request):
    scan = parse_float(request.form.get('scan', default=0), 0)
    scan = 0 if math.isnan(scan) else int(scan)
    thres = parse_float(request.form.get('thres', default=0.0), 0.0)
    mass = parse_float(request.form.get('mass', default=1.0), 1.0)
    clear = bool(request.form.get('clear', default=False))
    ext = request.form.get('ext', default='')
    predict = request.form.get('predict', default='{}')
    integration = request.form.get('integration', default='{}')
    multiplicity = request.form.get('multiplicity', default='{}')
    fname = parse_fname(request)
    simulatenrm = bool(request.form.get('simulatenrm', default=False))
    waveLength = request.form.get('wave_length', default=None)
    cyclicvolta = request.form.get('cyclic_volta', default=None)
    jcamp_idx = parse_int(request.form.get('jcamp_idx', default=0), 0)
    list_file_names = request.form.getlist('list_file_names[]')

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
        'integration': integration,
        'multiplicity': multiplicity,
        'fname': fname,
        'simulatenrm': simulatenrm,
        'waveLength': waveLength,
        'cyclic_volta': cyclicvolta,
        'jcamp_idx': jcamp_idx,
        'list_file_names': list_file_names,
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
        params.get('ext') or
        params.get('integration') or
        params.get('multiplicity') or
        params.get('fname') or
        params.get('simulatenrm')
    )
    if not has_params:
        params = False
    return params


def parse_array_to_dict_xys(peaks):
    xs, ys = peaks['x'], peaks['y']
    xys = []
    for idx, x in enumerate(xs):
        y = ys[idx]
        xys.append({'x': x, 'y': y})
    return xys
