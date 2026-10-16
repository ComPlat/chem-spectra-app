"""The band direction is decided once, from the quantity, and everything that
depends on it asks that one answer.

Two things still asked the technique instead, or forgot the answer between
requests:

- the auto-peak threshold had to be mapped through a conversion, and was --
  but only on the request that converted. A converted file coming back and
  re-picked (on `clear`, on a shift `delta`, or with no stored auto table)
  fell back to a linear mirror of the threshold, and lost its small peaks;
- the label headroom above the trace still keyed on the technique, so an
  infrared file in absorbance got the room a dip needs.
"""

import io
import zipfile

import pytest

from chem_spectra.lib.composer.technique import TechniqueComposer
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter


def _curve(tmp_path, datatype, yunits, ys, name, xunits='NANOMETERS'):
    xs = [400.0 + i for i in range(len(ys))]
    head = ['##TITLE=probe', '##JCAMP-DX=5.00', '##DATA TYPE=' + datatype,
            '##DATA CLASS=XYPOINTS', '##XUNITS=' + xunits,
            '##YUNITS=' + yunits,
            '##FIRSTX={}'.format(xs[0]), '##LASTX={}'.format(xs[-1]),
            '##MINX={}'.format(min(xs)), '##MAXX={}'.format(max(xs)),
            '##MINY={}'.format(min(ys)), '##MAXY={}'.format(max(ys)),
            '##NPOINTS={}'.format(len(xs)), '##FIRSTY={}'.format(ys[0]),
            '##XFACTOR=1', '##YFACTOR=1', '##XYPOINTS=(XY..XY)']
    rows = ['{:.1f}, {:.6f}'.format(x, y) for x, y in zip(xs, ys)]
    target = tmp_path / name
    target.write_text('\n'.join(head + rows + ['##END=']) + '\n')
    return target


def _chromatogram(tmp_path):
    """A 100 mAU band and a 12 mAU one on a 2 mAU baseline."""
    ys = [2.0] * 41
    ys[9], ys[10], ys[11] = 60.0, 100.0, 60.0
    ys[29], ys[30], ys[31] = 7.0, 12.0, 7.0
    return _curve(tmp_path, 'UV/VIS SPECTRUM', 'mAU', ys, 'chrom.jdx')


def _converter(path, **params):
    return JcampTechniqueConverter(JcampBaseConverter(str(path), params))


def _auto_x(converter):
    return sorted(round(float(x), 1) for x in converter.auto_peaks['x'])


@pytest.mark.parametrize('repick', [{'clear': True}, {'delta': 0.5}])
def test_a_converted_file_picks_the_same_peaks_when_it_comes_back(
        client, tmp_path, repick):
    """Converted once, then re-picked on the way back. Both bands, both
    times -- the second request has no `absorbance_range` to map through, so
    it has to recover the absorbance from the %T trace itself.
    """
    source = _chromatogram(tmp_path)
    converting = _converter(source, transmittance=True)
    assert _auto_x(converting) == [410.0, 430.0]

    with open(source, 'rb') as handle:
        response = client.post(
            '/zip_jcamp_n_img', content_type='multipart/form-data',
            data={'file': (io.BytesIO(handle.read()), 'chrom.jdx'),
                  'transmittance': 'true'})
    assert response.status_code == 200
    with zipfile.ZipFile(io.BytesIO(response.data)) as archive:
        name = next(n for n in archive.namelist() if n.endswith('.jdx'))
        composed = tmp_path / 'composed.jdx'
        composed.write_bytes(archive.read(name))

    returning = _converter(composed, **repick)
    assert returning.peaks_point_down is True
    assert returning.peak_threshold == pytest.approx(
        converting.peak_threshold, abs=1e-4)
    assert [round(x) for x in _auto_x(returning)] == [410, 430]


def _composer(tmp_path, datatype, yunits, ys, xunits='NANOMETERS', **params):
    path = _curve(tmp_path, datatype, yunits, ys, 'c.jdx', xunits=xunits)
    return TechniqueComposer(_converter(path, **params))


ABSORBANCE = [0.02] * 20 + [1.0]
PERCENT_T = [95.0] * 20 + [10.0]


@pytest.mark.parametrize('datatype, yunits, ys, xunits, params, headroom', [
    # the technique's own polarity: unchanged
    ('INFRARED SPECTRUM', '% TRANSMITTANCE', PERCENT_T, '1/CM', {}, 1.1),
    ('UV/VIS SPECTRUM', 'ABSORBANCE', ABSORBANCE, 'NANOMETERS', {}, 1.5),
    # the quantity reversing it
    ('INFRARED SPECTRUM', 'ABSORBANCE', ABSORBANCE, '1/CM', {}, 1.5),
    ('UV/VIS SPECTRUM', 'ABSORBANCE', ABSORBANCE, 'NANOMETERS',
     {'transmittance': True}, 1.1),
])
def test_the_label_headroom_follows_the_band_direction(
        tmp_path, datatype, yunits, ys, xunits, params, headroom):
    """Upward bands carry their labels above the trace and need the room;
    dips do not. Keyed on the technique, an infrared absorbance file got the
    dip's 1.1 and clipped its labels, and a converted UV/VIS trace the reverse.
    """
    composer = _composer(tmp_path, datatype, yunits, ys, xunits, **params)
    assert composer._TechniqueComposer__label_headroom() == headroom
