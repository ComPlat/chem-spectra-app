"""`/nmrium` converts a saved NMRium document to JCAMP — and is live.

The endpoint had no test and no fixture. It had both once: #259 ("remove unused
NMRium generation and harden parsing") deleted
`tests/lib/converter/nmrium/test_nmrium_base.py` and the fixtures
`nmrium_test_v3/v4/v5.nmrium`, leaving an empty `__pycache__` behind.

The premise was wrong. Traced on `chemotion_ELN@90c86acf2`:

    create_process (attachment_jcamp_aasm.rb:314, reached from :463)
      -> generate_spectrum_from_nmrium   when params[:ext] == 'nmrium'  (:317)
      -> CreateFromNMRium.jcamp_from_nmrium                             (:572)
      -> convert_nmrium_data -> POST {chemspectra}/nmrium      (jcamp.rb:484)

and `jcamp_files_already_present?` returns False for `.nmrium` (:155), so the
attachment is processed rather than skipped.

Not to be confused with the wrapper route: *viewing* a 2D spectrum hands the
zip or jcamp straight to nmrium-react-wrapper and this service is not involved.
*Saving* a `.nmrium` document is this path.

**The two document shapes both matter.** ELN PR #3586 wraps every save as
`{version, data}`; older documents are flat. `__read_file` branches on
`'data' in rawData`, and both fixtures here are real ELN documents — recovered
from `03e3fca^`, decimated to 1024 points so the file is 110 KB rather than
3 MB. Decimated, not truncated: the multiplicity ranges span the whole x axis,
and cutting the tail makes `__read_range_data` raise `IndexError` on an empty
slice.
"""

import json

import pytest
from werkzeug.datastructures import FileStorage

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.lib.converter.nmrium.base import NMRiumDataConverter

FLAT = './tests/fixtures/source/nmrium/flat.nmrium'
VERSIONED = './tests/fixtures/source/nmrium/versioned.nmrium'


def _convert(path):
    with open(path, 'rb') as handle:
        return NMRiumDataConverter(
            FileContainer(FileStorage(handle, filename='spectrum.nmrium')))


@pytest.fixture(params=[FLAT, VERSIONED], ids=['flat', 'versioned'])
def converted(request):
    return _convert(request.param)


# - - - the two document shapes - - -

def test_the_fixtures_really_are_the_two_shapes():
    """Guards the point of the pair. If both became wrapped, the flat case
    would silently stop being covered."""
    with open(FLAT) as handle:
        assert 'data' not in json.load(handle)
    with open(VERSIONED) as handle:
        assert 'data' in json.load(handle)


def test_both_shapes_yield_the_same_spectrum():
    """`{version, data}` is a wrapper, not a different document. ELN #3586
    makes every save wrapped, so the two must agree."""
    flat, versioned = _convert(FLAT), _convert(VERSIONED)
    assert flat.data['x'] == versioned.data['x']
    assert flat.data['y'] == versioned.data['y']


# - - - what the converter produces - - -

def test_the_curve_is_read(converted):
    assert converted.data is not None
    assert len(converted.data['x']) == len(converted.data['y']) == 1024


def test_it_declares_itself_nmr(converted):
    """No descriptor is carried, so since #297 this is the only converter
    relying on `BaseComposer._technique()`'s `non_nmr` fallback."""
    assert converted.non_nmr is False
    assert converted.is_2d is False


def test_datatype_comes_from_the_document_and_is_not_uppercased(converted):
    """Recorded because it is inconsistent, not because it is right.

    The constructor sets `datatype = 'NMR SPECTRUM'` and `datatypes =
    ['NMR SPECTRUM']`, then `__read_info` overwrites `datatype` from the
    document's `info.type` (`nmrium/base.py:158`) — NMRium writes
    `'NMR Spectrum'`. So the two attributes disagree in case on every upload,
    and nothing upper-cases this path the way `JcampBaseConverter.__set_datatype`
    canonicalises its own.

    Harmless today: `non_nmr` is False and the composer takes the NMR branch
    without consulting the string. It would stop being harmless the moment this
    core carried a descriptor looked up by datatype.
    """
    assert converted.datatype == 'NMR Spectrum'
    assert converted.datatypes == ['NMR SPECTRUM']


def test_the_x_axis_is_ppm_and_descending_in_range(converted):
    """A 1H spectrum: the boundary must span the real shift range, which is
    what decimating (rather than truncating) the fixture preserves."""
    assert converted.boundary['x']['min'] == pytest.approx(-4.01, abs=0.05)
    assert converted.boundary['x']['max'] == pytest.approx(15.99, abs=0.05)


def test_multiplicities_survive(converted):
    """The part that references x positions across the whole axis, and the
    reason the fixture is decimated rather than truncated: cutting the tail
    leaves a range with no points in it and `__read_range_data` raises
    `IndexError: index 0 is out of bounds` on the empty slice.
    """
    assert len(converted.mpy_itg_table) == 11
    for row in converted.mpy_itg_table:
        assert isinstance(row, str) and row.startswith('(')


def test_this_document_carries_no_edit_peaks(converted):
    """`edit_peaks` is a dict of parallel lists, not a list of peaks, so
    `len()` on it is 2 whatever it holds — which is how an earlier reading of
    this fixture mistook it for two peaks. It is empty.
    """
    assert converted.edit_peaks == {'x': [], 'y': []}


def test_a_datatable_is_produced(converted):
    """What the endpoint exists for: the ELN stores the JCAMP this becomes."""
    assert converted.datatable is not None
    assert len(converted.datatable) == 1024
    assert converted.datatable[0].strip().count(',') == 1, (
        'each row is "x, y" — this is what the composer writes into the JCAMP')


# - - - refusals - - -

def test_a_missing_file_is_not_a_crash():
    converter = NMRiumDataConverter(None)
    assert converter.data is None


def test_a_non_json_upload_is_not_a_crash(tmp_path):
    """`__read_file` catches ValueError and returns None, so the controller
    refuses rather than the converter raising."""
    junk = tmp_path / 'x.nmrium'
    junk.write_text('this is not json')
    assert _convert(str(junk)).data is None


# - - - through the endpoint, which is what the ELN actually calls - - -
#
# Everything above builds NMRiumDataConverter directly, so it would stay green
# if the route, `TraModel.tf_nmrium()`, the composer or the response handling
# broke. That is the gap this repository has been caught by twice: a converter
# test passing while the endpoint returned 500.

def _post(client, path):
    import io
    with open(path, 'rb') as handle:
        data = {'file': (io.BytesIO(handle.read()), 'spectrum.nmrium')}
    return client.post('/nmrium', content_type='multipart/form-data',
                       data=data)


@pytest.mark.parametrize('source', [FLAT, VERSIONED], ids=['flat', 'versioned'])
def test_the_endpoint_returns_a_jcamp_for_both_shapes(client, source):
    response = _post(client, source)
    assert response.status_code == 200
    body = response.data.decode('utf-8', 'ignore')
    assert body.startswith('##TITLE=')
    assert '##JCAMP-DX=' in body


@pytest.mark.parametrize('source', [FLAT, VERSIONED], ids=['flat', 'versioned'])
def test_the_endpoint_carries_the_curve_through(client, source):
    """The same 1024 points the converter reads, so a composer that silently
    dropped the data could not pass."""
    response = _post(client, source)
    body = response.data.decode('utf-8', 'ignore')
    assert '##NPOINTS=1024' in body
    assert '##XUNITS=ppm' in body
    assert '##FIRSTX=-4.0' in body and '##LASTX=15.9' in body
