"""A refused request leaves nothing behind for the next one.

Two ways it did. A refusal raised partway through an overlay left the curves
already drawn on the module-global pyplot figure, and the next image the same
worker rendered -- for a different request -- carried them. And an archive's
transmittance refusal came only after every other member had been converted
and rendered, which is work thrown away at best.
"""

import io
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

from chem_spectra.lib.converter.share import UnconvertibleSpectrum  # noqa: E402

MS_FILE = './tests/fixtures/source/ms/svs813f1.jdx'


def _curve(datatype, yunits='ABSORBANCE', xunits='NANOMETERS'):
    """A small valid spectrum; absorbance, so it converts when asked."""
    ys = [0.02] * 20 + [1.0]
    xs = [200.0 + i for i in range(len(ys))]
    head = ['##TITLE=probe', '##JCAMP-DX=5.00', '##DATA TYPE=' + datatype,
            '##DATA CLASS=XYPOINTS', '##XUNITS=' + xunits,
            '##YUNITS=' + yunits,
            '##FIRSTX={}'.format(xs[0]), '##LASTX={}'.format(xs[-1]),
            '##MINX={}'.format(min(xs)), '##MAXX={}'.format(max(xs)),
            '##MINY={}'.format(min(ys)), '##MAXY={}'.format(max(ys)),
            '##NPOINTS={}'.format(len(xs)), '##FIRSTY={}'.format(ys[0]),
            '##XFACTOR=1', '##YFACTOR=1', '##XYPOINTS=(XY..XY)']
    rows = ['{:.1f}, {:.4f}'.format(x, y) for x, y in zip(xs, ys)]
    return '\n'.join(head + rows + ['##END=']) + '\n'


def test_a_refused_overlay_leaves_no_curve_for_the_next_request(client):
    """The UV/VIS file sorts first and is drawn; the mass spectrum then
    refuses the transmittance request. Before, that curve stayed on the
    shared figure and turned up in whatever this worker drew next.
    """
    plt.close('all')
    with open(MS_FILE, 'rb') as handle:
        ms = handle.read()
    response = client.post(
        '/combine_images', content_type='multipart/form-data',
        data={'files[]': [(io.BytesIO(_curve('UV/VIS SPECTRUM').encode()),
                           'a_uv.jdx'),
                          (io.BytesIO(ms), 'b_ms.jdx')],
              'transmittance': 'true'})
    assert response.status_code == 422
    assert 'does not measure absorption' in json.loads(response.data)['error']
    assert len(plt.gca().lines) == 0, (
        'a refused overlay left %d curve(s) on the shared figure'
        % len(plt.gca().lines))


def test_an_archive_refuses_transmittance_before_converting_anything(
        tmp_path, monkeypatch):
    """The LC/MS group is known as soon as the members have been read, so the
    refusal can come then -- not after every other member was converted and
    rendered.
    """
    import chem_spectra.lib.converter.bagit.base as bagit_base

    bag = tmp_path / 'bag'
    data = bag / 'data'
    data.mkdir(parents=True)
    # sorts first, and is no part of the LC/MS group
    (data / 'a_ir.jdx').write_text(
        _curve('INFRARED SPECTRUM', xunits='1/CM'))
    # with the chromatogram, these two make the archive an LC/MS dataset
    (data / 'b_uv.jdx').write_text(_curve('UV/VIS SPECTRUM'))
    (data / 'c_trace.jdx').write_text(_curve('HPLC UV/VIS SPECTRUM'))
    (bag / 'bagit.txt').write_text('BagIt-Version: 0.97\n')

    composed = []

    def spy(*args, **kwargs):
        composed.append(args)
        raise AssertionError('a member was composed before the refusal')

    monkeypatch.setattr(bagit_base, 'TechniqueComposer', spy)
    with pytest.raises(UnconvertibleSpectrum, match='LC/MS dataset'):
        bagit_base.BagItBaseConverter(str(bag), params={'transmittance': True})
    assert composed == []
