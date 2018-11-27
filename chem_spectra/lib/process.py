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
