"""Which Bruker uploads are refused as multi-dimensional, and what /file/convert
answers for a processed 1D upload.

The layouts are synthetic: the refusal is decided from file names before
anything is read, so empty files are enough.
"""
import base64
import io
import os
import zipfile

import pytest

from chem_spectra.lib.converter.share import UnconvertibleSpectrum
from chem_spectra.model.transformer import refuse_if_bruker_multi_dimensional


def _tree(root, paths):
    for path in paths:
        full = os.path.join(root, path)
        os.makedirs(os.path.dirname(full), exist_ok=True)
        open(full, 'wb').close()


def _zip(paths):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w') as z:
        for path in paths:
            z.writestr(path, b'')
    return buf.getvalue()


# a DOSY (2D, `ser`) carrying its 1D reference proton experiment inside it:
# the layout of #250's upload
DOSY_WITH_NESTED_1D = [
    '151/acqus', '151/acqu2s', '151/ser',
    '151/150/acqus', '151/150/fid',
]


def test_a_1d_experiment_inside_a_2d_one_is_refused(tmp_path):
    _tree(tmp_path, DOSY_WITH_NESTED_1D)
    with pytest.raises(UnconvertibleSpectrum) as err:
        refuse_if_bruker_multi_dimensional(str(tmp_path))
    assert str(err.value) == (
        '151: this is a 2D NMR experiment. ChemSpectra reads one-dimensional '
        'spectra only; open it in NMRium instead')


def test_a_2d_experiment_alone_is_refused_without_a_prefix(tmp_path):
    _tree(tmp_path, ['acqus', 'acqu2s', 'ser'])
    with pytest.raises(UnconvertibleSpectrum, match='^this is a 2D NMR'):
        refuse_if_bruker_multi_dimensional(str(tmp_path))


def test_a_3d_experiment_says_3d(tmp_path):
    _tree(tmp_path, ['7/acqus', '7/acqu2s', '7/acqu3s', '7/ser'])
    with pytest.raises(UnconvertibleSpectrum, match='^7: this is a 3D NMR'):
        refuse_if_bruker_multi_dimensional(str(tmp_path))


def test_a_1d_experiment_beside_a_2d_one_is_not_refused(tmp_path):
    # a sample folder with a proton spectrum next to an HSQC: the 1D one is
    # what gets converted, and the 2D one does not stop it
    _tree(tmp_path, ['sample/10/acqus', 'sample/10/fid',
                     'sample/12/acqus', 'sample/12/acqu2s', 'sample/12/ser'])
    refuse_if_bruker_multi_dimensional(str(tmp_path))


def test_a_1d_upload_is_not_refused(tmp_path):
    _tree(tmp_path, ['1/acqus', '1/fid', '1/pdata/1/1r'])
    refuse_if_bruker_multi_dimensional(str(tmp_path))


@pytest.mark.parametrize('url', [
    '/zip_jcamp_n_img', '/api/v1/chemspectra/file/convert'])
def test_the_endpoints_refuse_it_with_the_reason(client, url):
    response = client.post(
        url, data={'file': (io.BytesIO(_zip(DOSY_WITH_NESTED_1D)), '151.zip')},
        content_type='multipart/form-data')
    assert response.status_code == 422
    assert response.json['error'].startswith('151: this is a 2D NMR experiment')


def test_file_convert_returns_every_spectrum_of_a_processed_bruker_upload(client):
    # the transform of the raw FID and Bruker's processed spectrum, as
    # /zip_jcamp_n_img returns them; this used to be an HTML 500
    with open('./tests/fixtures/source/bruker/1H.zip', 'rb') as f:
        response = client.post(
            '/api/v1/chemspectra/file/convert',
            data={'file': (f, '1H.zip')}, content_type='multipart/form-data')
    assert response.status_code == 200
    jcamps = response.json['list_jcamps']
    assert len(jcamps) == 2
    for jcamp in jcamps:
        text = base64.b64decode(jcamp).decode('utf-8', errors='ignore')
        assert '##DATA TYPE=NMR SPECTRUM' in text or '##DATA TYPE=LINK' in text
