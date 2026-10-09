"""A block with no readable x range is refused, not failed.

`JcampTechniqueConverter.__read_xs` builds the x axis from FIRSTX/LASTX, or
from FIRST/LAST over the observe frequency. When neither gives both ends it
built the axis from `None + delta`: a TypeError and a 500. That happened to an
NTUPLES block read by `##PAGE=` (whichever of UV-VIS or HPLC UV-VIS it was
labelled) and to any file without a range. It is now an
`UnconvertibleSpectrum` whose reason says which of the two it is: a 422 on the
conversion endpoints, and the outline on `/predict/infrared`, which reports its
errors in the body.
"""

import io
import json

import pytest

IR = './tests/fixtures/source/IR.dx'
MOLFILE = './tests/fixtures/source/molfile/svs813f1_B.mol'

PAGED_NTUPLES = """##TITLE=Spectrum
##JCAMP-DX=5.00 $$ chemotion-converter-app (1.8.0)
##DATA TYPE={datatype}
##DATA CLASS=NTUPLES
##XUNITS=MINUTES
##YUNITS=SIGNAL
##NPOINTS=2
##PAGE=Wavelength= 210.0
##XYDATA=(XY..XY)
0.0, 1.0
1.0, 2.0
##END=
"""


def _read(path):
    with open(path) as handle:
        return handle.read()


def _ir_without_range():
    """IR.dx with its ##FIRSTX/##LASTX records removed."""
    return '\n'.join(line for line in _read(IR).splitlines()
                     if not line.startswith(('##FIRSTX', '##LASTX'))) + '\n'


def _post(client, route, text, name='paged.jdx'):
    data = {'file': (io.BytesIO(text.encode()), name)}
    return client.post(route, content_type='multipart/form-data', data=data)


@pytest.mark.parametrize('datatype', ['UV-VIS', 'HPLC UV-VIS'])
@pytest.mark.parametrize('route', [
    '/api/v1/chemspectra/file/convert',
    '/zip_jcamp_n_img',
])
def test_a_paged_block_is_refused_with_a_reason(client, route, datatype):
    response = _post(client, route, PAGED_NTUPLES.format(datatype=datatype))
    assert response.status_code == 422
    assert 'paged NTUPLES' in json.loads(response.data)['error']


@pytest.mark.parametrize('route', [
    '/api/v1/chemspectra/file/convert',
    '/zip_jcamp_n_img',
])
def test_a_file_without_a_range_is_refused_without_naming_ntuples(client, route):
    """An ordinary spectrum missing its range is not told it is a paged block."""
    response = _post(client, route, _ir_without_range(), 'ir.dx')
    assert response.status_code == 422
    error = json.loads(response.data)['error']
    assert 'x range' in error
    assert 'NTUPLES' not in error


def test_a_half_reading_does_not_hide_a_valid_firstx_lastx(client):
    """FIRST without LAST used to set one end and skip the FIRSTX/LASTX
    reading, so a file with a valid range was refused. Each reading now gives
    both ends or nothing."""
    text = _read(IR).replace(
        '##XUNITS=', '##.OBSERVE FREQUENCY=100.0\n##FIRST=3997.453, 0\n##XUNITS=', 1)
    assert '##FIRST=' in text and '##LAST=' not in text
    response = _post(client, '/zip_jcamp_n_img', text, 'ir_first_only.dx')
    assert response.status_code == 200


def test_predict_infrared_reports_it_in_the_outline(client):
    """This endpoint reports errors in the body: the refusal must not turn its
    200 into a 422 (the ELN reads only a 200), and the outline carries the
    reason."""
    response = client.post(
        '/predict/infrared',
        data={
            'spectrum': (io.BytesIO(_ir_without_range().encode()), 'x.jdx'),
            'molfile': (io.BytesIO(_read(MOLFILE).encode()), 'm.mol'),
        },
        content_type='multipart/form-data',
    )
    assert response.status_code == 200
    outline = response.get_json()['outline']
    assert outline['code'] == 400
    assert 'x range' in outline['text']


def test_combine_images_refuses_with_the_reason(client):
    response = client.post(
        '/combine_images',
        data={'files[]': [
            (io.BytesIO(_read(IR).encode()), 'a.dx'),
            (io.BytesIO(PAGED_NTUPLES.format(datatype='UV-VIS').encode()), 'b.jdx'),
        ]},
        content_type='multipart/form-data',
    )
    assert response.status_code == 422
    assert 'paged NTUPLES' in response.get_json()['error']


def test_a_block_with_an_x_range_is_untouched(client):
    response = _post(client, '/zip_jcamp_n_img', _read(IR), 'IR.dx')
    assert response.status_code == 200
