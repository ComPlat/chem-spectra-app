import io
import zipfile
from os.path import basename


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
            fname = get_fname(abs_path, filename)
            zf.write(abs_path, fname)
    memory.seek(0)
    for tmp in tmp_arr:
        tmp.close()
    return memory


def extract_params(request):
    params = {
        'peaks_str': request.form.get('peaks_str', default=None),
        'select_x': request.form.get('shift_select_x', default=None),
        'ref_name': request.form.get('shift_ref_name', default=None),
        'ref_value': request.form.get('shift_ref_value', default=None),
    }
    has_params = (
        params.get('peaks_str') or
        params.get('select_x') or
        params.get('ref_name') or
        params.get('ref_value')
    )
    if not has_params:
        params = False
    return params
