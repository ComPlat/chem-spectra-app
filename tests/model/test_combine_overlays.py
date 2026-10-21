"""A combined image annotates one spectrum: the one `jcamp_idx` names.

Peaks, integrals, multiplets and cyclic voltammetry max/min pairs are drawn
by `TechniqueComposer.plot_overlays`, for that spectrum alone and without
narrowing the x range the other curves share. These tests spy on that call
(per-file images call it too, with `adjust_xlim=True`, so only the overlay's
calls are kept) and on the axes at the moment the image is saved.
"""

import io
import json
import zipfile

import pytest
from werkzeug.datastructures import FileStorage

import chem_spectra.lib.converter.bagit.base as bagit_module
import chem_spectra.model.transformer as transformer_module
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.lib.composer.technique import TechniqueComposer
from chem_spectra.lib.converter.bagit.base import BagItBaseConverter
from chem_spectra.model.transformer import TransformerModel

H1 = './tests/fixtures/source/1H.dx'
C13 = './tests/fixtures/source/13C-CPD.dx'
# one block, so a record added under its first line is the target block's:
# 1H.dx and 13C-CPD.dx are LINK files, where it would land in the root
SINGLE_BLOCK_NMR = './tests/fixtures/source/SVS-790A_13C.jdx'
CHROMATOGRAM = './tests/fixtures/source/hplc/chromatogram.jdx'
CV = './tests/fixtures/source/cyclicvoltammetry/RCV_LSH-R444_full+Fc.jdx'
CV_BAGIT = './tests/fixtures/source/bagit/cv/File053_BagIt.zip'


def _read(path):
    with open(path, 'rb') as handle:
        return handle.read()


def _inverted(body):
    first, rest = body.split(b'\n', 1)
    return first + b'\n##$CSINVERTY=true\n' + rest


@pytest.fixture
def overlays(monkeypatch):
    """Every overlay drawn onto a combined image, as the composer saw it."""
    calls = []
    real = TechniqueComposer.plot_overlays

    def spy(self, plt, y_values=None, adjust_xlim=True):
        if not adjust_xlim:
            calls.append({
                'core': self.core,
                'ncl': self.core.ncl,
                'draw_y_inverted': getattr(self.core, 'draw_y_inverted', False),
                'list_max_min_peaks': self.core.params.get('list_max_min_peaks'),
                'cv_density_scale': self._cv_density_scale,
            })
        return real(self, plt, y_values, adjust_xlim)

    monkeypatch.setattr(TechniqueComposer, 'plot_overlays', spy)
    return calls


@pytest.fixture
def saved_axes(monkeypatch):
    captured = {}
    real_savefig = transformer_module.plt.savefig

    def spy(*args, **kwargs):
        ax = transformer_module.plt.gca()
        captured['xlim'], captured['ylim'] = ax.get_xlim(), ax.get_ylim()
        return real_savefig(*args, **kwargs)

    monkeypatch.setattr(transformer_module.plt, 'savefig', spy)
    return captured


def _combine(bodies, names, params=None, extras=None):
    files = [
        FileContainer(FileStorage(io.BytesIO(body), filename=name))
        for body, name in zip(bodies, names)
    ]
    model = TransformerModel(None, params={'ext': 'dx', **(params or {})},
                             multiple_files=files)
    return model.tf_combine(list_file_names=names, extraParams=extras)


def test_the_overlay_goes_on_the_spectrum_sent_at_jcamp_idx(overlays):
    """`jcamp_idx` counts files as sent; the combine sorts them by name."""
    tf = _combine([_read(H1), _read(C13)], ['b.dx', 'a.dx'], {'jcamp_idx': 0})
    tf.close()
    assert [call['ncl'] for call in overlays] == ['1H']


def test_files_sharing_a_name_get_one_overlay(overlays):
    tf = _combine([_read(H1), _read(C13)], ['same.dx', 'same.dx'],
                  {'jcamp_idx': 1})
    tf.close()
    assert [call['ncl'] for call in overlays] == ['13C']


def test_a_string_index_is_tolerated(overlays):
    # form params arrive untyped when the controller's parsing is bypassed
    tf = _combine([_read(H1), _read(C13)], ['a.dx', 'b.dx'], {'jcamp_idx': '1'})
    tf.seek(0)
    assert tf.read(8) == b'\x89PNG\r\n\x1a\n'
    tf.close()
    assert [call['ncl'] for call in overlays] == ['13C']


@pytest.mark.parametrize('jcamp_idx', ['5', -1, 'x'])
def test_an_index_naming_no_file_still_combines(overlays, jcamp_idx):
    tf = _combine([_read(H1), _read(H1)], ['a.dx', 'b.dx'],
                  {'jcamp_idx': jcamp_idx})
    tf.seek(0)
    assert tf.read(8) == b'\x89PNG\r\n\x1a\n'
    tf.close()
    # 'x' parses as 0; the others name no file and draw no overlay
    assert len(overlays) == (1 if jcamp_idx == 'x' else 0)


def test_an_nmr_overlay_reads_high_to_low(overlays, saved_axes):
    _combine([_read(H1), _read(H1)], ['a.dx', 'b.dx']).close()
    left, right = saved_axes['xlim']
    assert left > right


def test_an_overlay_with_a_forward_curve_reads_low_to_high(overlays, saved_axes):
    _combine([_read(H1), _read(CHROMATOGRAM)], ['a.dx', 'b.jdx']).close()
    left, right = saved_axes['xlim']
    assert left < right


def test_the_overlay_does_not_narrow_the_shared_x_range(overlays, saved_axes):
    """Peak labels on a 1H spectrum pin its view to -0.2..9 ppm; the 13C
    curve beside it runs to 200."""
    _combine([_read(H1), _read(C13)], ['a.dx', 'b.dx'],
             {'jcamp_idx': 0, 'peaks_str': '1.1,1.112#2.2,2.224',
              'ref_name': 'DMSO-d6 (quin)', 'ref_value': '2.5'}).close()
    assert len(overlays) == 1
    assert max(saved_axes['xlim']) > 150


def test_a_mixed_overlay_draws_its_inverted_active_spectrum_upright(
        overlays, saved_axes):
    """The axis stays upright when the curves disagree, so the active
    spectrum's markers and labels must not point as if it were flipped."""
    nmr = _read(SINGLE_BLOCK_NMR)
    _combine([_inverted(nmr), nmr], ['a.jdx', 'b.jdx'], {'jcamp_idx': 0}).close()
    low, high = saved_axes['ylim']
    assert low < high
    assert [call['draw_y_inverted'] for call in overlays] == [False]
    # and the converter is left as the file declared it
    assert overlays[0]['core'].draw_y_inverted is True


def test_an_inverted_overlay_draws_its_active_spectrum_inverted(
        overlays, saved_axes):
    inverted = _inverted(_read(SINGLE_BLOCK_NMR))
    _combine([inverted, inverted], ['a.jdx', 'b.jdx'], {'jcamp_idx': 1}).close()
    top, bottom = saved_axes['ylim']
    assert top > bottom
    assert [call['draw_y_inverted'] for call in overlays] == [True]


CV_PEAKS = [
    {'max': {'x': 0.5, 'y': 1e-5}, 'min': {'x': 0.3, 'y': -1e-5},
     'isRef': True, 'pecker': {'x': 0.4, 'y': 5e-6}},
    {'max': {'x': 1.0, 'y': 2e-5}, 'min': {'x': 0.8, 'y': -2e-5},
     'isRef': False},
]


def test_cv_peaks_sent_in_extras_reach_the_active_spectrum(overlays):
    """The combine endpoint may send the CV state in `extras` alone."""
    state = {
        'spectraList': [{'list': [], 'shift': {}, 'hasRefPeak': False},
                        {'list': CV_PEAKS, 'shift': {}, 'hasRefPeak': True}],
        'cvDisplay': {'yScaleFactor': 2},
    }
    tf = _combine([_read(CV), _read(CV)], ['a.jdx', 'b.jdx'], {'jcamp_idx': 1},
                  extras=json.dumps({'cyclicvolta': json.dumps(state)}))
    tf.close()
    assert len(overlays) == 1
    assert overlays[0]['list_max_min_peaks'] == CV_PEAKS
    assert overlays[0]['cv_density_scale'] == 2.0


def _bagit(tmp_path, params, invert=()):
    with zipfile.ZipFile(CV_BAGIT) as archive:
        archive.extractall(tmp_path)
    for name in invert:
        member = tmp_path / 'data' / name
        member.write_bytes(_inverted(member.read_bytes()))
    return BagItBaseConverter(str(tmp_path), params)


def test_the_bagit_overlay_goes_on_the_member_at_jcamp_idx(tmp_path, overlays):
    converter = _bagit(tmp_path, {'jcamp_idx': 1})
    assert converter.combined_image is not None
    assert len(overlays) == 1
    assert overlays[0]['core'] is converter._composers[1].core


def test_a_mixed_bagit_overlay_restores_the_member_it_drew(
        tmp_path, overlays, monkeypatch):
    captured = {}
    real_savefig = bagit_module.plt.savefig

    def spy(*args, **kwargs):
        captured['ylim'] = bagit_module.plt.gca().get_ylim()
        return real_savefig(*args, **kwargs)

    monkeypatch.setattr(bagit_module.plt, 'savefig', spy)
    converter = _bagit(tmp_path, {'jcamp_idx': 0}, invert=('table_01.jdx',))
    low, high = captured['ylim']
    assert low < high
    assert [call['draw_y_inverted'] for call in overlays] == [False]
    # the composer is kept on the converter; it must not keep the override
    assert converter._composers[0].core.draw_y_inverted is True


# - - - the ELN's save and refresh paths - - -
#
# Both convert each file of `dst_list` with `jcamp_idx` set to that file's
# own index, so each picks its own CV state and axis units. The combined
# image is built after that loop, and used to inherit the last index: the
# overlay landed on the last spectrum, carrying the peaks the user had edited
# on the one they were viewing.

@pytest.mark.parametrize('route', [
    '/api/v1/chemspectra/file/save',
    '/api/v1/chemspectra/file/refresh',
])
@pytest.mark.parametrize('jcamp_idx, ncl', [('0', '1H'), ('1', '13C')])
def test_save_and_refresh_annotate_the_requested_spectrum(
        client, overlays, route, jcamp_idx, ncl):
    data = {
        'src': (io.BytesIO(_read(H1)), 'a.dx'),
        'dst_list': [(io.BytesIO(_read(H1)), 'a.dx'),
                     (io.BytesIO(_read(C13)), 'b.dx')],
        'jcamp_idx': jcamp_idx,
    }
    response = client.post(route, data=data, content_type='multipart/form-data')
    assert response.status_code == 200
    assert [call['ncl'] for call in overlays] == [ncl]
