"""A JCAMP nmrglue cannot extract a data array from.

`make_ni_data_ys` dereferenced `base.data.shape` unguarded. nmrglue returns
None when it parses no data array, so such a file raised AttributeError
straight out of the request -- a 500 for what is really just an unusable
upload.

It now raises UnparsableJcampData, which jcamp2cvp turns into the same
falsy result any other unconvertible file produces, so the controllers'
existing `abort(...)` handles it. No new contract was invented: this is the
status those endpoints already return when conversion yields nothing.
"""

import io

import pytest

from chem_spectra import create_app
from chem_spectra.lib.converter.jcamp.data_parse import UnparsableJcampData

# 129 data-block markers, but nmrglue yields data=None
source_unparsable = './tests/fixtures/source/MS.dx'


@pytest.fixture
def client():
    return create_app({'TESTING': True}).test_client()


def _post(client, endpoint):
    with open(source_unparsable, 'rb') as handle:
        return client.post(
            endpoint,
            data={'file': (io.BytesIO(handle.read()), 'unparsable.jdx')},
            content_type='multipart/form-data',
        )


@pytest.mark.parametrize('endpoint,expected', [
    # each endpoint's own existing "could not convert" status, unchanged
    ('/api/v1/chemspectra/file/convert', 400),
    ('/zip_jcamp_n_img', 403),
])
def test_unparsable_jcamp_is_rejected_not_a_500(client, endpoint, expected):
    """Asserted at the endpoint. A converter-level test would pass while
    the controller still returned 500 -- that has happened twice."""
    assert _post(client, endpoint).status_code == expected


def test_the_parse_failure_is_named(client):
    """The exception says what went wrong, rather than AttributeError on
    'NoneType' from three frames down."""
    from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
    from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter

    base = JcampBaseConverter(source_unparsable)
    assert base.data is None
    with pytest.raises(UnparsableJcampData):
        JcampNIConverter(base)


def test_a_parsable_file_is_unaffected(client):
    """The guard must not catch working files."""
    with open('./tests/fixtures/source/IR.dx', 'rb') as handle:
        response = client.post(
            '/api/v1/chemspectra/file/convert',
            data={'file': (io.BytesIO(handle.read()), 'IR.dx')},
            content_type='multipart/form-data',
        )
    assert response.status_code == 200


# - - - the other controller paths that build a converter - - -
#
# Guarding jcamp2cvp alone was not enough: five other sites construct a
# JcampNIConverter, and four of them are reachable from an upload. Found by
# review of #294, after the first version of this change claimed the class
# of bug was handled when only one path was.

source_molfile = './tests/fixtures/source/molfile/svs813f1_B.mol'
source_good = './tests/fixtures/source/IR.dx'


def _bad():
    return io.BytesIO(open(source_unparsable, 'rb').read())


def _good():
    return io.BytesIO(open(source_good, 'rb').read())


def test_predict_by_peaks_form_rejects_it(client):
    """to_converter() now returns False, and the caller dereferenced it.

    `cv.edit_peaks` on a bool raised AttributeError
    (controller/inference_api.py).
    """
    response = client.post(
        '/predict/by_peaks_form',
        data={
            'spectrum': (_bad(), 'x.jdx'),
            'molfile': (io.BytesIO(open(source_molfile, 'rb').read()), 'm.mol'),
            'layout': '13C',
            # 'peaks' omitted on purpose: the handler defaults it to '{}'.
            # Passing '' makes json.loads blow up before the code under
            # test is reached -- a separate robustness gap, not this one.
        },
        content_type='multipart/form-data',
    )
    assert response.status_code == 400


def test_predict_infrared_reports_it_in_the_outline(client):
    """This endpoint reports errors in the body, not the status code."""
    response = client.post(
        '/predict/infrared',
        data={
            'spectrum': (_bad(), 'x.jdx'),
            'molfile': (io.BytesIO(open(source_molfile, 'rb').read()), 'm.mol'),
        },
        content_type='multipart/form-data',
    )
    assert response.status_code == 200
    assert response.get_json()['outline']['code'] == 400


def test_combine_images_skips_an_unusable_file(client):
    """One bad file must not lose the overlay of the good ones."""
    response = client.post(
        '/combine_images',
        data={'files[]': [(_bad(), 'a.jdx'), (_good(), 'b.dx')]},
        content_type='multipart/form-data',
    )
    assert response.status_code == 200


def test_combine_images_rejects_an_overlay_of_nothing(client):
    """All files unusable: an empty image would look like success."""
    response = client.post(
        '/combine_images',
        data={'files[]': [(_bad(), 'a.jdx'), (_bad(), 'b.jdx')]},
        content_type='multipart/form-data',
    )
    assert response.status_code == 400


def test_combine_images_is_unaffected_for_usable_files(client):
    response = client.post(
        '/combine_images',
        data={'files[]': [(_good(), 'a.dx'), (_good(), 'b.dx')]},
        content_type='multipart/form-data',
    )
    assert response.status_code == 200
