"""The y-axis units we write back, and what it takes to change them.

ChemSpectra used to decide two things about a spectrum's y-axis without
consulting each other or the person who sent the file:

- whether to invert the data -- from its *shape* (`median < 0.5 * max`);
- what units to write -- from the *declared* string containing "absorb".

So an infrared file declaring `A.U.` had its numbers mirrored and its label
left alone, and an HPLC UV/VIS file declaring `ABSORBANCE` came back saying
`TRANSMITTANCE` -- the inverse quantity -- with the untouched original still
visible as `###YUNITS`. The file contradicted itself.

Both decisions now belong to the client, as an explicit instruction, and
nothing is inferred: `transmittance` converts, `T = 100 * 10**(-A)`, and says
so. Without it, `##YUNITS` is exactly what arrived.

`invert_y` is the second such instruction, and it is a **drawing** one. It
used to mirror the stored array with `max(y) - y` -- a display concern
implemented as a data mutation, which destroys the baseline, skews the
area-under-curve integration, travels into the exported JCAMP, and had to be
labelled `X - inverted`, a unit JCAMP-DX does not define. It now moves the
viewport instead: `plt.ylim` is drawn the other way up, the numbers are
untouched, and the request is recorded as `##$CSINVERTY` for renderers
downstream.

The need is regional rather than exotic -- DSC's exo-up against exo-down,
cyclic voltammetry's IUPAC against Texas sign convention, and indirect
photometric HPLC, where the analyte is a dip. Instrument software offers the
same toggle. What none of them do is rewrite the stored trace.
"""

import io
import json
import re

import numpy as np
import pytest

import chem_spectra.lib.composer.technique as technique_module
from chem_spectra.lib.composer.technique import TechniqueComposer
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import (
    JcampTechniqueConverter, UnconvertibleSpectrum, is_transmittance_unit,
)
from chem_spectra.lib.converter.jcamp.techniques import SPECTRUM_TECHNIQUES

ABSORBANCE_SHAPED = './tests/fixtures/source/JPK-948.jdx'   # declares ABSORBANCE
TRANSMITTANCE_SHAPED = './tests/fixtures/source/IR.dx'      # declares TRANSMITTANCE


def _probe(source, datatype, tmp_path, yunits=None, params=None):
    body = open(source).read()
    body = body.replace(
        re.search(r'##DATA TYPE=.*', body).group(0), '##DATA TYPE=' + datatype, 1)
    if yunits is not None:
        body = body.replace(
            re.search(r'##YUNITS=.*', body).group(0), '##YUNITS=' + yunits, 1)
    target = tmp_path / 'probe.jdx'
    target.write_text(body)
    return JcampTechniqueConverter(JcampBaseConverter(str(target), params))


def _absorbance_probe(tmp_path, datatype='INFRARED SPECTRUM', params=None):
    """A spectrum whose y really is absorbance: a low baseline with bands up to
    A = 2. No fixture has this -- JPK-948 declares ABSORBANCE but carries
    detector counts in the thousands, which the conversion rightly refuses.
    """
    xs = [400.0 + i for i in range(101)]
    ys = [0.02] * 101
    for centre, height in ((20, 2.0), (55, 1.2), (80, 0.6)):
        for offset, scale in ((-1, 0.4), (0, 1.0), (1, 0.4)):
            ys[centre + offset] = height * scale
    body = [
        '##TITLE=synthetic absorbance\n',
        '##JCAMP-DX=5.00\n',
        '##DATA TYPE=' + datatype + '\n',
        '##DATA CLASS=XYPOINTS\n',
        '##FIRSTX={}\n'.format(xs[0]),
        '##LASTX={}\n'.format(xs[-1]),
        '##MINX={}\n'.format(min(xs)),
        '##MAXX={}\n'.format(max(xs)),
        '##MINY={}\n'.format(min(ys)),
        '##MAXY={}\n'.format(max(ys)),
        '##NPOINTS={}\n'.format(len(xs)),
        '##FIRSTY={}\n'.format(ys[0]),
        '##XUNITS=1/CM\n',
        '##YUNITS=ABSORBANCE\n',
        '##XYPOINTS=(XY..XY)\n',
    ]
    body += ['{:.6f}, {:.6f}\n'.format(x, y) for x, y in zip(xs, ys)]
    body.append('##END=\n')
    target = tmp_path / 'absorbance.jdx'
    target.write_text(''.join(body))
    return JcampTechniqueConverter(JcampBaseConverter(str(target), params))


# - - - with no instructions, nothing is touched - - -

@pytest.mark.parametrize('key', sorted(SPECTRUM_TECHNIQUES))
def test_declared_units_are_returned_unchanged(key, tmp_path):
    """Every technique in the registry, so one added later is covered too."""
    if key == 'MS':
        pytest.skip('MS has its own converter and composer')
    converter = _probe(ABSORBANCE_SHAPED, key, tmp_path, yunits='ABSORBANCE')
    assert converter.label['y'] == 'ABSORBANCE'


@pytest.mark.parametrize('declared', ['ABSORBANCE', 'A.U.', 'ARBITRARY', '%T'])
def test_infrared_is_not_special_cased_any_more(declared, tmp_path):
    """`ABSORBANCE` and `A.U.` were the two that went wrong, in opposite ways."""
    converter = _probe(TRANSMITTANCE_SHAPED, 'INFRARED SPECTRUM', tmp_path,
                       yunits=declared)
    assert converter.label['y'] == declared
    assert list(converter.ys) == list(converter.data)


def test_hplc_uv_vis_keeps_its_absorbance(tmp_path):
    """The case that prompted this: HPLC UV/VIS is always absorbance."""
    converter = _probe(ABSORBANCE_SHAPED, 'HPLC UV/VIS SPECTRUM', tmp_path,
                       yunits='ABSORBANCE')
    assert converter.label['y'] == 'ABSORBANCE'


def test_composed_units_agree_with_the_preserved_original(tmp_path):
    """The defect's real consequence, and what nothing asserted before."""
    converter = _probe(ABSORBANCE_SHAPED, 'HPLC UV/VIS SPECTRUM', tmp_path,
                       yunits='ABSORBANCE')
    meta = ''.join(TechniqueComposer(converter).meta)
    ours = re.search(r'^##YUNITS=(.*)$', meta, re.M).group(1).strip()
    original = re.search(r'^###YUNITS=(.*)$', meta, re.M).group(1).strip()
    assert ours == original == 'ABSORBANCE'


# - - - invert: the viewport, never the array - - -

def _rendered_ylim(converter):
    """The axis limits as drawn. They have to be read during `savefig`:
    `tf_img` clears the figure straight after, which is why
    `test_ms_orientation.py` spies in the same place."""
    captured = {}
    real_savefig = technique_module.plt.savefig

    def spy(*args, **kwargs):
        captured['ylim'] = technique_module.plt.gca().get_ylim()
        return real_savefig(*args, **kwargs)

    technique_module.plt.savefig = spy
    try:
        TechniqueComposer(converter).tf_img().close()
    finally:
        technique_module.plt.savefig = real_savefig
    return captured['ylim']


def test_invert_leaves_the_numbers_alone(tmp_path):
    """The whole point. Same array, same units."""
    plain = _absorbance_probe(tmp_path)
    inverted = _absorbance_probe(tmp_path, params={'invert_y': True})
    assert list(inverted.ys) == list(plain.ys)
    assert inverted.label['y'] == plain.label['y']


def test_invert_does_not_invent_a_unit(tmp_path):
    """`ABSORBANCE - inverted` was a label JCAMP-DX cannot parse."""
    inverted = _absorbance_probe(tmp_path, params={'invert_y': True})
    assert 'inverted' not in inverted.label['y']


def test_the_mutation_is_gone_from_the_source():
    """Not merely unreachable from the API."""
    source = open('chem_spectra/lib/converter/jcamp/technique.py').read()
    assert 'np.max(ys) - ys' not in source
    assert '- inverted' not in source


def test_invert_draws_the_axis_the_other_way_up(tmp_path):
    """Measured from the render rather than asserted off the flag."""
    low, high = _rendered_ylim(_absorbance_probe(tmp_path))
    assert low < high, 'the plain render should ascend'

    top, bottom = _rendered_ylim(
        _absorbance_probe(tmp_path, params={'invert_y': True}))
    assert top > bottom, 'the inverted render should descend'
    assert (top, bottom) == pytest.approx((high, low)), (
        'only the direction should change, not the bounds')


def test_invert_is_recorded_for_renderers_downstream(tmp_path):
    """`##$CSINVERTY` is a viewing preference now, not a claim about the
    data -- which is what makes it safe to write."""
    composed = TechniqueComposer(
        _absorbance_probe(tmp_path, params={'invert_y': True}))
    assert '##$CSINVERTY=true' in ''.join(composed.meta)


def test_invert_and_transmittance_are_no_longer_exclusive(tmp_path):
    """They were, while inversion meant `max(%T) - %T` -- fractional
    absorptance, not linear in concentration and unlabellable. Nothing
    computes that any more: one converts the numbers, the other turns the
    picture over, so asking for both is coherent."""
    converter = _absorbance_probe(
        tmp_path, params={'transmittance': True, 'invert_y': True})
    assert converter.label['y'] == '% TRANSMITTANCE'
    assert min(converter.ys) == pytest.approx(1.0, abs=1e-4)


# - - - transmittance - - -

def test_transmittance_converts_and_relabels(tmp_path):
    """Percent, not the 0-1 ratio.

    %T is how instruments commonly present it; a 0-1 array is routinely
    misread as absorbance, which itself runs 0-2.5. `% TRANSMITTANCE` is not a
    JCAMP-DX unit -- 4.24 defines TRANSMITTANCE only as the ratio I_T/I_0 --
    so this is a readability choice, not spec conformance.
    """
    converter = _absorbance_probe(tmp_path, params={'transmittance': True})
    assert converter.label['y'] == '% TRANSMITTANCE'
    # T% = 100 * 10**(-A): the A = 2.0 band becomes 1.0, the 0.02 baseline ~95.5
    assert min(converter.ys) == pytest.approx(1.0, abs=1e-4)
    assert max(converter.ys) == pytest.approx(100 * 10 ** -0.02, abs=1e-4)


def test_transmittance_refuses_a_file_that_declares_transmittance(tmp_path):
    """What the file says outranks what its shape suggests.

    `IR.dx` declares `##YUNITS=TRANSMITTANCE`, so this is refused on a fact
    rather than on the median heuristic below -- which matters because that
    heuristic misfires on a heavily absorbing sample.
    """
    with pytest.raises(UnconvertibleSpectrum, match='already declares'):
        _probe(TRANSMITTANCE_SHAPED, 'INFRARED SPECTRUM', tmp_path,
               params={'transmittance': True})


@pytest.mark.parametrize('yunits', [
    'TRANSMITTANCE', '% TRANSMITTANCE', '%T', 'T', 'T%', 'TRANSMISSION',
    'Transmission (%)', 'transmittance [%]', 'Percent Transmittance',
])
def test_every_common_spelling_of_transmittance_is_refused(yunits, tmp_path):
    """Instrument exports rarely say exactly 'TRANSMITTANCE'.

    The data is absorbance-shaped on purpose -- median far below the maximum,
    as for a heavily absorbing %T trace -- so the median heuristic would let
    every one of these through, and only the unit can stop the conversion.
    """
    ys = [0.1] * 100
    ys[40:60] = [2.0] * 20
    with pytest.raises(UnconvertibleSpectrum, match='already declares'):
        JcampTechniqueConverter(JcampBaseConverter(
            _synthetic(tmp_path, ys, yunits=yunits), {'transmittance': True}))


@pytest.mark.parametrize('yunits', [
    'ABSORBANCE', 'Absorbance (a.u.)', 'ARBITRARY UNITS', 'TEMPERATURE',
    'OPTICAL DENSITY',
])
def test_units_that_are_not_transmittance_still_convert(yunits, tmp_path):
    """An exact set, so no unrelated unit blocks a conversion by accident."""
    ys = [0.1] * 100
    ys[40:60] = [2.0] * 20
    converter = JcampTechniqueConverter(JcampBaseConverter(
        _synthetic(tmp_path, ys, yunits=yunits), {'transmittance': True}))
    assert converter.label['y'] == '% TRANSMITTANCE'


def test_a_single_units_record_is_read_for_the_declared_unit(tmp_path):
    """JCAMP 6 declares x, y and z in one `##UNITS=` record per block.

    Review caught this: the guard read `UNITS[1]`, copied from `__set_label`,
    where index 1 is the spectrum block of an NMR LINK file. Four fixtures
    (the MNOVA and MS v6 files) declare exactly one record and no `##YUNITS=`
    at all, so a hardcoded 1 found nothing and a file declaring transmittance
    only there was converted a second time. The target block's index is the
    right one, and for a single-block file that is 0.
    """
    xs = [400.0 + i for i in range(60)]
    ys = [95.0] * 60
    ys[20] = 2.0
    body = [
        '##TITLE=single units record\n', '##JCAMP-DX=6.00\n',
        '##DATA TYPE=INFRARED SPECTRUM\n', '##DATA CLASS=XYPOINTS\n',
        '##UNITS=1/CM, % TRANSMITTANCE, ARBITRARY\n',
        '##FIRSTX={}\n'.format(xs[0]), '##LASTX={}\n'.format(xs[-1]),
        '##MINX={}\n'.format(min(xs)), '##MAXX={}\n'.format(max(xs)),
        '##MINY={}\n'.format(min(ys)), '##MAXY={}\n'.format(max(ys)),
        '##NPOINTS={}\n'.format(len(xs)), '##FIRSTY={}\n'.format(ys[0]),
        '##XYPOINTS=(XY..XY)\n',
    ]
    body += ['{:.6f}, {:.6f}\n'.format(x, y) for x, y in zip(xs, ys)]
    body.append('##END=\n')
    target = tmp_path / 'single_units.jdx'
    target.write_text(''.join(body))

    with pytest.raises(UnconvertibleSpectrum, match='already declares'):
        JcampTechniqueConverter(
            JcampBaseConverter(str(target), {'transmittance': True}))


def _single_units_record(tmp_path, yunits, units_y):
    """A one-block file declaring its y unit twice: ##YUNITS= and a JCAMP 6
    ##UNITS= triple, which overrides it."""
    ys = [0.1] * 100
    ys[40:60] = [2.0] * 20
    path = _synthetic(tmp_path, ys, yunits=yunits)
    body = open(path).read().replace(
        '##XYPOINTS=', '##UNITS=1/CM, {}, ARBITRARY UNITS\n##XYPOINTS='.format(
            units_y), 1)
    open(path, 'w').write(body)
    return path


def test_the_guard_and_the_label_read_the_same_record(tmp_path):
    """Review caught this: the guard read UNITS[target_idx] then UNITS[0],
    __set_label read UNITS[1]. In a one-record file the guard refused on the
    ##UNITS= transmittance while the composed label came from ##YUNITS=."""
    path = _single_units_record(tmp_path, 'ABSORBANCE', '% TRANSMITTANCE')
    label = JcampTechniqueConverter(JcampBaseConverter(path, None)).label
    assert 'TRANSMITTANCE' in label['y'].upper()
    with pytest.raises(UnconvertibleSpectrum, match='already declares'):
        JcampTechniqueConverter(
            JcampBaseConverter(path, {'transmittance': True}))


def test_the_units_record_overrides_yunits_for_both(tmp_path):
    """The other direction: ##UNITS= says absorbance, ##YUNITS= says %T.
    Both readers take the triple, so the label and the guard agree that this
    is absorbance and the conversion goes ahead."""
    path = _single_units_record(tmp_path, '% TRANSMITTANCE', 'ABSORBANCE')
    assert JcampTechniqueConverter(
        JcampBaseConverter(path, None)).label['y'] == 'ABSORBANCE'
    converted = JcampTechniqueConverter(
        JcampBaseConverter(path, {'transmittance': True}))
    assert converted.label['y'] == '% TRANSMITTANCE'


def test_the_shape_guard_still_catches_a_file_that_declares_nothing(tmp_path):
    """The heuristic remains, as a second line for files that declare a unit
    naming neither quantity. Transmittance-shaped: baseline near the maximum.
    """
    ys = [0.95] * 100
    ys[20] = 0.01
    with pytest.raises(UnconvertibleSpectrum, match='already appears'):
        JcampTechniqueConverter(JcampBaseConverter(
            _synthetic(tmp_path, ys, yunits='ARBITRARY'),
            {'transmittance': True}))


def test_transmittance_refuses_data_that_is_not_absorbance_scaled(tmp_path):
    """JPK-948 declares ABSORBANCE but carries detector counts of 13-6168.
    10**(-6168) underflows, so converting would return a flat line labelled
    TRANSMITTANCE -- worse than anything this change fixes."""
    with pytest.raises(UnconvertibleSpectrum, match='not absorbance'):
        _probe(ABSORBANCE_SHAPED, 'INFRARED SPECTRUM', tmp_path,
               yunits='ABSORBANCE',
               params={'transmittance': True, 'jcamp_idx': 0})


# - - - provenance - - -

def test_the_flags_are_emitted_only_when_asked_for(tmp_path):
    plain = TechniqueComposer(_absorbance_probe(tmp_path))
    assert '##$CSTRANSMITTANCE' not in ''.join(plain.meta)
    assert '##$CSINVERTY' not in ''.join(plain.meta)

    converted = TechniqueComposer(
        _absorbance_probe(tmp_path, params={'transmittance': True}))
    assert '##$CSTRANSMITTANCE=true' in ''.join(converted.meta)


def test_nothing_infers_from_the_data_shape():
    """The 0.5 heuristic must not come back as a decision."""
    source = open('chem_spectra/lib/converter/jcamp/technique.py').read()
    read_ys = source[source.index('def __read_ys'):source.index('def __to_transmittance')]
    assert '0.5' not in read_ys


# - - - at the endpoint, which is where it matters - - -

def _post(client, source, **form):
    import io
    with open(source, 'rb') as handle:
        data = dict(file=(io.BytesIO(handle.read()), source.split('/')[-1]), **form)
    return client.post('/zip_jcamp_n_img', content_type='multipart/form-data',
                       data=data)


def test_endpoint_without_flags_succeeds(client):
    assert _post(client, TRANSMITTANCE_SHAPED).status_code == 200


def test_endpoint_invert_succeeds(client):
    """Through the controller, because a converter-level test would pass
    while the controller still returned 500; that has happened twice here."""
    assert _post(client, TRANSMITTANCE_SHAPED, invert_y='true').status_code == 200


def test_endpoint_accepts_both_instructions_together(client):
    """`IR.dx` declares transmittance, so `transmittance` is still refused on
    that fact. The reason names the declared unit rather than a conflict
    between the two instructions, because there is no longer a conflict."""
    response = _post(client, TRANSMITTANCE_SHAPED,
                     transmittance='true', invert_y='true')
    assert response.status_code == 422
    assert 'already declares' in json.loads(response.data)['error']


@pytest.mark.parametrize('source, reason', [
    (TRANSMITTANCE_SHAPED, 'already declares'),
    (ABSORBANCE_SHAPED, 'not absorbance'),
])
def test_endpoint_refuses_an_impossible_conversion(client, source, reason):
    """422 with the reason, not a 500 and not a ruined spectrum.

    A converter-level test would pass while the controller still returned 500 --
    that has happened twice in this repo.
    """
    response = _post(client, source, transmittance='true')
    assert response.status_code == 422
    assert reason in response.get_json()['error']


# - - - the conversion must never emit non-finite values (#298 review) - - -

def _synthetic(tmp_path, ys, yunits='ABSORBANCE'):
    xs = [400.0 + i for i in range(len(ys))]
    body = [
        '##TITLE=synthetic\n', '##JCAMP-DX=5.00\n',
        '##DATA TYPE=INFRARED SPECTRUM\n', '##DATA CLASS=XYPOINTS\n',
        '##FIRSTX={}\n'.format(xs[0]), '##LASTX={}\n'.format(xs[-1]),
        '##MINX={}\n'.format(min(xs)), '##MAXX={}\n'.format(max(xs)),
        '##MINY={}\n'.format(min(ys)), '##MAXY={}\n'.format(max(ys)),
        '##NPOINTS={}\n'.format(len(xs)), '##FIRSTY={}\n'.format(ys[0]),
        '##XUNITS=1/CM\n', '##YUNITS={}\n'.format(yunits),
        '##XYPOINTS=(XY..XY)\n',
    ]
    body += ['{:.6f}, {:.6f}\n'.format(x, y) for x, y in zip(xs, ys)]
    body.append('##END=\n')
    target = tmp_path / 'synthetic.jdx'
    target.write_text(''.join(body))
    return str(target)


def test_large_negative_values_are_refused_not_turned_into_infinity(tmp_path):
    """`10**(400)` overflows to inf, and inf passed both original guards.

    A trace at -400 clears the median check (its median is not near its max)
    and the ceiling check (its max is negative), so the request used to
    succeed with non-finite data and bounds.
    """
    ys = [-400.0] * 100
    ys[:5] = [-1.0] * 5
    path = _synthetic(tmp_path, ys)
    with pytest.raises(UnconvertibleSpectrum, match='overflow'):
        JcampTechniqueConverter(
            JcampBaseConverter(path, {'transmittance': True}))


def test_slightly_negative_absorbance_still_converts(tmp_path):
    """Baseline drift below zero is normal and must not be refused."""
    ys = [-0.05] * 100
    ys[20] = 1.5
    converter = JcampTechniqueConverter(
        JcampBaseConverter(_synthetic(tmp_path, ys), {'transmittance': True}))
    assert converter.label['y'] == '% TRANSMITTANCE'
    assert all(abs(float(v)) < 1e6 for v in converter.ys)


# - - - the records have to survive a round trip - - -
#
# Every pass through this app recomposes, and a recompose carries no
# instructions, so a record written on one request was dropped on the next.
# That made both flags write-only: the editor, which reads them from the file,
# never saw one. Found by loading a composed file in the editor and watching
# the y-axis toggle stay off.

def _recompose(client, body):
    response = client.post(
        '/jcamp',
        data={'file': (io.BytesIO(body), 'p.jdx')},
        content_type='multipart/form-data',
    )
    assert response.status_code == 200
    return response.data


def _records(body, name):
    return re.findall((r'^##\$%s=.*$' % name).encode(), body, re.M)


def _recompose_with(client, body, **form):
    data = {'file': (io.BytesIO(body), 'p.jdx')}
    data.update(form)
    response = client.post('/jcamp', data=data,
                           content_type='multipart/form-data')
    assert response.status_code == 200
    return response.data


def test_the_inversion_record_survives_recomposition(client):
    with open(TRANSMITTANCE_SHAPED, 'rb') as handle:
        source = handle.read()
    once = _recompose_with(client, source, invert_y='true')
    assert _records(once, 'CSINVERTY'), 'the first pass must write it'
    twice = _recompose(client, once)
    assert _records(twice, 'CSINVERTY'), 'a recompose must keep it'
    assert _records(_recompose(client, twice), 'CSINVERTY')


def test_a_file_that_never_asked_does_not_gain_the_record(client):
    """The other half: reading the file must not invent a flag."""
    with open(TRANSMITTANCE_SHAPED, 'rb') as handle:
        plain = _recompose(client, handle.read())
    assert not _records(plain, 'CSINVERTY')


@pytest.mark.parametrize('record, form', [
    ('CSINVERTY', {'invert_y': 'true'}),
    ('CSTRANSMITTANCE', {'transmittance': 'true'}),
])
def test_the_records_are_not_echoed_into_the_metadata(
        client, tmp_path, record, form):
    """They are written as real records when they apply. Echoed as well,
    `###$CSINVERTY= true` survived the record being cleared, and the file's
    metadata contradicted itself."""
    _absorbance_probe(tmp_path)
    once = _recompose_with(
        client, (tmp_path / 'absorbance.jdx').read_bytes(), **form)
    twice = _recompose(client, once)
    assert _records(twice, record)
    assert ('###$' + record).encode() not in twice


def test_clearing_the_inversion_leaves_no_trace_of_it(client):
    with open(TRANSMITTANCE_SHAPED, 'rb') as handle:
        once = _recompose_with(client, handle.read(), invert_y='true')
    cleared = _recompose_with(client, _recompose(client, once),
                              invert_y='false')
    assert b'CSINVERTY' not in cleared


def test_an_old_echo_neither_sets_the_record_nor_survives(client):
    """`###$CSINVERTY= true` is metadata, not a record: nmrglue keys it as
    `#$CSINVERTY`. Files composed before the echo was suppressed carry it, so
    it must not set the flag, and is dropped on the next pass."""
    with open(TRANSMITTANCE_SHAPED, 'rb') as handle:
        source = handle.read()
    first, rest = source.split(b'\n', 1)
    legacy = first + b'\n###$CSINVERTY= true\n' + rest
    recomposed = _recompose(client, legacy)
    assert b'CSINVERTY' not in recomposed


def test_the_transmittance_record_survives_recomposition(client, tmp_path):
    """Through the endpoint, which is the path that actually recomposes.

    `_absorbance_probe` is called only to write the file; the conversion under
    test is the one the request asks for.
    """
    _absorbance_probe(tmp_path)
    source = (tmp_path / 'absorbance.jdx').read_bytes()
    once = _recompose_with(client, source, transmittance='true')
    assert _records(once, 'CSTRANSMITTANCE')
    assert b'##YUNITS=% TRANSMITTANCE' in once
    twice = _recompose(client, once)
    assert _records(twice, 'CSTRANSMITTANCE'), 'a recompose must keep it'
    assert b'##YUNITS=% TRANSMITTANCE' in twice


def test_an_explicit_false_clears_the_inversion_record(client):
    """The record is a preference the request may override in either
    direction. `invert_y` used to collapse "not sent" and "false" into one
    value, so once a file was inverted nothing could un-invert it."""
    with open(TRANSMITTANCE_SHAPED, 'rb') as handle:
        once = _recompose_with(client, handle.read(), invert_y='true')
    assert _records(once, 'CSINVERTY')
    cleared = _recompose_with(client, once, invert_y='false')
    assert not _records(cleared, 'CSINVERTY')
    assert not _records(_recompose(client, cleared), 'CSINVERTY')


@pytest.mark.parametrize('value', ['undefined', 'null', 'on'])
def test_an_unrecognised_value_keeps_the_inversion_record(client, value):
    """Only a recognised false clears it. `undefined` is what JS FormData
    makes of a missing value, and must not undo the user's choice."""
    with open(TRANSMITTANCE_SHAPED, 'rb') as handle:
        once = _recompose_with(client, handle.read(), invert_y='true')
    kept = _recompose_with(client, once, invert_y=value)
    assert _records(kept, 'CSINVERTY')


def _converted_once(client, tmp_path):
    _absorbance_probe(tmp_path)
    source = (tmp_path / 'absorbance.jdx').read_bytes()
    return _recompose_with(client, source, transmittance='true')


def _absorbance_units():
    return json.dumps({'axes': [{'xUnit': '', 'yUnit': 'ABSORBANCE'}]})


def test_a_recompose_cannot_relabel_converted_data(client, tmp_path):
    """The record outlives the request that converted, and so must the label:
    otherwise %T data goes out declared as absorbance while still carrying
    ##$CSTRANSMITTANCE=true."""
    once = _converted_once(client, tmp_path)
    twice = _recompose_with(client, once, axes_units=_absorbance_units())
    assert _records(twice, 'CSTRANSMITTANCE')
    assert b'##YUNITS=% TRANSMITTANCE' in twice
    assert b'##YUNITS=ABSORBANCE' not in twice


def test_converted_data_is_refused_a_second_conversion(client, tmp_path):
    """The record says the data is already transmittance, so that is the
    reason given -- not whatever the shape heuristic happens to conclude."""
    once = _converted_once(client, tmp_path)
    response = client.post(
        '/zip_jcamp_n_img', content_type='multipart/form-data',
        data={'file': (io.BytesIO(once), 'p.jdx'), 'transmittance': 'true'})
    assert response.status_code == 422
    assert 'already converted' in response.get_json()['error']


def test_a_relabelled_conversion_is_still_refused(client, tmp_path):
    """The record, not the unit text, is what proves the data is %T.

    A converted file whose y unit was later rewritten (here to `COUNTS`,
    which names no transmittance) slips past the declared-unit check. Without
    the record check a strongly absorbing trace would be converted a second
    time, and anything else refused with a misleading shape reason.
    """
    once = _converted_once(client, tmp_path)
    relabelled = re.sub(rb'^##YUNITS=.*$', b'##YUNITS=COUNTS', once,
                        flags=re.M)
    assert _records(relabelled, 'CSTRANSMITTANCE')
    # every declared y unit must be one the declared-unit check lets through
    declared = re.findall(rb'^##YUNITS=(.*)$', relabelled, re.M)
    assert declared
    assert not any(is_transmittance_unit(u.decode()) for u in declared)
    assert not re.search(rb'^##UNITS=', relabelled, re.M)
    response = client.post(
        '/zip_jcamp_n_img', content_type='multipart/form-data',
        data={'file': (io.BytesIO(relabelled), 'p.jdx'),
              'transmittance': 'true'})
    assert response.status_code == 422
    assert 'already converted' in response.get_json()['error']


# - - - the overlay image draws the same way up as the single one - - -
#
# Under #298 the mirrored data reached tf_combine by itself. Now inversion is
# a viewport flip, and the overlay has its own figure, so it has to apply the
# flip too -- otherwise an inverted file is drawn inverted alone and upright
# beside others.

def _with_inversion_record(source):
    first, rest = source.split(b'\n', 1)
    return first + b'\n##$CSINVERTY=true\n' + rest


def _overlay_ylim(bodies):
    import chem_spectra.model.transformer as transformer_module
    from werkzeug.datastructures import FileStorage
    from chem_spectra.controller.helper.file_container import FileContainer
    from chem_spectra.model.transformer import TransformerModel

    files = [
        FileContainer(FileStorage(io.BytesIO(body), filename='s%d.dx' % idx))
        for idx, body in enumerate(bodies)
    ]
    captured = {}
    real_savefig = transformer_module.plt.savefig

    def spy(*args, **kwargs):
        captured['ylim'] = transformer_module.plt.gca().get_ylim()
        return real_savefig(*args, **kwargs)

    transformer_module.plt.savefig = spy
    try:
        TransformerModel(None, params={'ext': 'dx'},
                         multiple_files=files).tf_combine().close()
    finally:
        transformer_module.plt.savefig = real_savefig
    return captured['ylim']


def test_the_overlay_flips_when_every_curve_is_inverted():
    with open(TRANSMITTANCE_SHAPED, 'rb') as handle:
        plain = handle.read()
    inverted = _with_inversion_record(plain)
    low, high = _overlay_ylim([plain, plain])
    assert low < high, 'the plain overlay should ascend'
    top, bottom = _overlay_ylim([inverted, inverted])
    assert (top, bottom) == pytest.approx((high, low)), (
        'only the direction should change, not the bounds')


def test_a_mixed_overlay_stays_upright():
    """One y axis per figure: flipping it would misdraw the upright curves."""
    with open(TRANSMITTANCE_SHAPED, 'rb') as handle:
        plain = handle.read()
    low, high = _overlay_ylim([plain, _with_inversion_record(plain)])
    assert low < high


def _bagit_overlay_ylim(tmp_path, invert):
    """`invert` names which of the archive's members carry the record."""
    import zipfile
    import chem_spectra.lib.converter.bagit.base as bagit_module
    from chem_spectra.lib.converter.bagit.base import BagItBaseConverter

    with zipfile.ZipFile(
            './tests/fixtures/source/bagit/cv/File053_BagIt.zip') as archive:
        archive.extractall(tmp_path)
    for name in invert:
        member = tmp_path / 'data' / name
        member.write_bytes(_with_inversion_record(member.read_bytes()))

    captured = {}
    real_savefig = bagit_module.plt.savefig

    def spy(*args, **kwargs):
        captured['ylim'] = bagit_module.plt.gca().get_ylim()
        return real_savefig(*args, **kwargs)

    bagit_module.plt.savefig = spy
    try:
        BagItBaseConverter(str(tmp_path))
    finally:
        bagit_module.plt.savefig = real_savefig
    return captured['ylim']


BAGIT_MEMBERS = ('table_01.jdx', 'table_02.jdx', 'table_03.jdx')


def test_the_bagit_overlay_flips_when_every_member_is_inverted(tmp_path):
    """The archive overlay is a separate figure from tf_combine's."""
    low, high = _bagit_overlay_ylim(tmp_path / 'plain', ())
    assert low < high, 'the plain overlay should ascend'
    top, bottom = _bagit_overlay_ylim(tmp_path / 'inverted', BAGIT_MEMBERS)
    assert (top, bottom) == pytest.approx((high, low))


def test_a_mixed_bagit_overlay_stays_upright(tmp_path):
    low, high = _bagit_overlay_ylim(tmp_path, BAGIT_MEMBERS[:1])
    assert low < high


# - - - unit records belong to the block that declares them - - -
#
# nmrglue collects every block's records into one list per label, so an
# index into it is a block number only when every block declares the record.
# Here only the interferogram declares ##UNITS=, and its `CM, VOLTS, ...`
# landed at index 0 -- the spectrum's index -- relabelling the spectrum and
# hiding a %T y axis from the transmittance guard.

LINK_WITH_INTERFEROGRAM = './tests/fixtures/source/ir_link_interferogram.jdx'


def test_another_blocks_units_do_not_relabel_the_spectrum():
    converter = JcampTechniqueConverter(
        JcampBaseConverter(LINK_WITH_INTERFEROGRAM, None))
    assert converter.label == {'x': '1/CM', 'y': 'ABSORBANCE'}


def test_the_composed_file_keeps_the_spectrums_units():
    composed = ''.join(TechniqueComposer(JcampTechniqueConverter(
        JcampBaseConverter(LINK_WITH_INTERFEROGRAM, None))).meta)
    assert '##XUNITS=1/CM' in composed
    assert '##YUNITS=ABSORBANCE' in composed
    assert 'VOLTS' not in composed.split('$$ === CHEMSPECTRA SPECTRUM ORIG')[0]


def test_the_guard_reads_the_spectrums_own_unit(tmp_path):
    """The data is absorbance-shaped, so only the unit can stop this."""
    body = open(LINK_WITH_INTERFEROGRAM).read().replace(
        '##YUNITS=ABSORBANCE', '##YUNITS=%T', 1)
    target = tmp_path / 'link.jdx'
    target.write_text(body)
    with pytest.raises(UnconvertibleSpectrum, match='already declares'):
        JcampTechniqueConverter(
            JcampBaseConverter(str(target), {'transmittance': True}))


def test_a_trailing_comma_in_the_units_record_is_accepted(tmp_path):
    """Mnova writes `HZ, ARBITRARY UNITS, ARBITRARY UNITS,`. Rejected, the
    label was taken from the peak-table block's XUNITS/YUNITS instead."""
    converter = JcampTechniqueConverter(JcampBaseConverter(
        './tests/fixtures/source/mnova/STM212_H.jcamp', None))
    assert converter.label == {'x': 'PPM', 'y': 'ARBITRARY'}


# - - - peak tables are converted with the trace - - -
#
# The normal ELN flow saves a file once and asks for the conversion later, so
# the file already carries peak tables picked on the absorbance trace. Kept
# as read, PEAK TABLE AUTO stayed in absorbance (y = 0.02 against a 1-95 %T
# trace) and the preview drew its markers along the bottom.

def _composed_absorbance(client, tmp_path, **form):
    _absorbance_probe(tmp_path)
    return _recompose_with(
        client, (tmp_path / 'absorbance.jdx').read_bytes(), **form)


def _read_back(tmp_path, body):
    target = tmp_path / 'read_back.jdx'
    target.write_bytes(body)
    return JcampTechniqueConverter(JcampBaseConverter(str(target), None))


def test_auto_peaks_are_picked_again_on_the_converted_trace(client, tmp_path):
    saved = _composed_absorbance(client, tmp_path)
    converted = _read_back(
        tmp_path, _recompose_with(client, saved, transmittance='true'))
    low, high = float(min(converted.ys)), float(max(converted.ys))
    assert high > 90, 'the trace should be %T'
    # the absorbance bands at 420, 455 and 480 are the %T dips
    assert sorted(converted.auto_peaks['x']) == [420.0, 455.0, 480.0]
    assert all(low <= y <= high for y in converted.auto_peaks['y'])
    assert min(converted.auto_peaks['y']) == pytest.approx(1.0)


def test_edited_peaks_keep_their_position_and_change_units(client, tmp_path):
    saved = _composed_absorbance(client, tmp_path, peaks_str='420,2.0#455,1.2')
    converted = _read_back(
        tmp_path, _recompose_with(client, saved, transmittance='true'))
    assert converted.edit_peaks['x'] == [420.0, 455.0]
    assert converted.edit_peaks['y'] == pytest.approx(
        [100 * 10 ** -2.0, 100 * 10 ** -1.2])


def test_peaks_sent_with_the_conversion_are_converted_too(tmp_path):
    """The editor sending them was showing the absorbance trace."""
    converter = _absorbance_probe(
        tmp_path, params={'transmittance': True, 'peaks_str': '480,0.6'})
    assert converter.edit_peaks['x'] == [480.0]
    assert converter.edit_peaks['y'] == pytest.approx([100 * 10 ** -0.6])


def test_peaks_are_untouched_without_a_conversion(tmp_path):
    converter = _absorbance_probe(tmp_path, params={'peaks_str': '480,0.6'})
    assert converter.edit_peaks['y'] == [0.6]


def test_integrals_in_the_file_refuse_the_conversion(tmp_path):
    _absorbance_probe(tmp_path)
    source = tmp_path / 'absorbance.jdx'
    source.write_text(source.read_text().replace(
        '##XYPOINTS=',
        '##$OBSERVEDINTEGRALS= (X Y Z)\n(425.0, 415.0, 1.0)\n##XYPOINTS=', 1))
    with pytest.raises(UnconvertibleSpectrum, match='integrals'):
        JcampTechniqueConverter(
            JcampBaseConverter(str(source), {'transmittance': True}))


def test_integrals_in_the_request_refuse_the_conversion(tmp_path):
    stack = json.dumps({'stack': [{'xL': 415.0, 'xU': 425.0, 'area': 1.0}],
                        'refArea': 1, 'refFactor': 1, 'shift': 0})
    with pytest.raises(UnconvertibleSpectrum, match='integrals'):
        _absorbance_probe(
            tmp_path, params={'transmittance': True, 'integration': stack})


def test_an_empty_integral_table_does_not_refuse(tmp_path):
    """Only rows count: the record's first line is its `(X Y Z)` header."""
    _absorbance_probe(tmp_path)
    source = tmp_path / 'absorbance.jdx'
    source.write_text(source.read_text().replace(
        '##XYPOINTS=', '##$OBSERVEDINTEGRALS= (X Y Z)\n##XYPOINTS=', 1))
    converter = JcampTechniqueConverter(
        JcampBaseConverter(str(source), {'transmittance': True}))
    assert converter.converted_to_transmittance


MULTI_BLOCK = 'tests/fixtures/source/ir_link_interferogram.jdx'


def _composed_jcamp(response):
    """The .jdx out of the endpoint's zip."""
    import zipfile
    with zipfile.ZipFile(io.BytesIO(response.data)) as archive:
        name = next(n for n in archive.namelist() if n.endswith('.jdx'))
        return archive.read(name).decode('utf-8', errors='ignore')


def test_the_endpoint_labels_a_multi_block_file_from_its_own_block(client):
    """The unit records come from the block that holds the spectrum.

    At the endpoint the uploaded file is a NamedTemporaryFile that is closed,
    and so deleted, before the technique converter is built. A reader that
    re-opens the path by name therefore found nothing here and fell back to
    nmrglue's flattened lists, where the interferogram's `CM`/`VOLTS` sit
    beside the spectrum's own records. The same file was then labelled one way
    through a fixture path and another way through the endpoint.
    """
    response = _post(client, MULTI_BLOCK)
    assert response.status_code == 200
    composed = _composed_jcamp(response)
    assert '##XUNITS=1/CM' in composed
    assert '##YUNITS=ABSORBANCE' in composed


def test_a_byte_order_mark_does_not_shift_the_blocks(client, tmp_path):
    """`str.strip` keeps U+FEFF, so a BOM'd first line does not start with
    `##`. The outer `##TITLE=` was then missed, every block moved up one, and
    the spectrum was labelled from the interferogram that follows it."""
    source = tmp_path / 'bom.jdx'
    source.write_text(open(MULTI_BLOCK, encoding='utf-8').read(),
                      encoding='utf-8-sig')
    composed = _composed_jcamp(_post(client, str(source)))
    assert '##XUNITS=1/CM' in composed
    assert '##YUNITS=ABSORBANCE' in composed
