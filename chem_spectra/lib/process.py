import io
import zipfile


def to_zip_response(tmp_arr):
    memory = io.BytesIO()
    with zipfile.ZipFile(memory, 'w') as zf:
        for tmp in tmp_arr:
            zf.write(tmp.name)
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
