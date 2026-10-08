"""A block with no readable x axis is refused, not failed.

An NTUPLES block read by `##PAGE=` declares no FIRSTX/LASTX (nor FIRST/LAST),
so `JcampTechniqueConverter.__read_xs` found no range and built the axis from
`None + delta`: a TypeError and a 500, on master, whichever of UV-VIS or
HPLC UV-VIS the block was labelled. It is now an `UnconvertibleSpectrum`,
which every endpoint answers with a 422 and the reason.
"""

import io
import json

import pytest

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


def _post(client, route, text, name='paged.jdx'):
    data = {'file': (io.BytesIO(text.encode()), name)}
    return client.post(route, content_type='multipart/form-data', data=data)


@pytest.mark.parametrize('datatype', ['UV-VIS', 'HPLC UV-VIS'])
@pytest.mark.parametrize('route', [
    '/api/v1/chemspectra/file/convert',
    '/zip_jcamp_n_img',
])
def test_a_block_without_an_x_axis_is_refused_with_a_reason(client, route,
                                                            datatype):
    response = _post(client, route, PAGED_NTUPLES.format(datatype=datatype))
    assert response.status_code == 422
    assert 'no x axis' in json.loads(response.data)['error']


def test_a_block_with_an_x_range_is_untouched(client):
    response = _post(client, '/zip_jcamp_n_img',
                     open('./tests/fixtures/source/IR.dx').read(), 'IR.dx')
    assert response.status_code == 200
