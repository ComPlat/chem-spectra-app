"""The y-axis units we write back, and what it takes to change them.

ChemSpectra used to decide two things about a spectrum's y-axis without
consulting each other or the person who sent the file:

- whether to invert the data -- from its *shape* (`median < 0.5 * max`);
- what units to write -- from the *declared* string containing "absorb".

So an infrared file declaring `A.U.` had its numbers mirrored and its label
left alone, and an HPLC UV/VIS file declaring `ABSORBANCE` came back saying
`TRANSMITTANCE` -- the inverse quantity -- with the untouched original still
visible as `###YUNITS`. The file contradicted itself.

Both decisions now belong to the client, as explicit instructions, and
nothing is inferred:

- `transmittance` converts, `T = 10**(-A)`, and says so;
- `invert` mirrors for display and says so.

With neither, `##YUNITS` is exactly what arrived.
"""

import re

import pytest

from chem_spectra.lib.composer.technique import TechniqueComposer
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import (
    JcampTechniqueConverter, UnconvertibleSpectrum,
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


# - - - invert - - -

def test_invert_mirrors_and_says_so(tmp_path):
    converter = _probe(ABSORBANCE_SHAPED, 'INFRARED SPECTRUM', tmp_path,
                       yunits='ABSORBANCE', params={'invert_y': True})
    assert list(converter.ys) == [max(converter.data) - y for y in converter.data]
    assert converter.label['y'] == 'ABSORBANCE - inverted'


def test_invert_is_not_gated_to_infrared(tmp_path):
    """It is an instruction; refusing it for other techniques would be silent."""
    converter = _probe(ABSORBANCE_SHAPED, 'HPLC UV/VIS SPECTRUM', tmp_path,
                       yunits='mAU', params={'invert_y': True})
    assert converter.label['y'] == 'mAU - inverted'


# - - - transmittance - - -

def test_transmittance_converts_and_relabels(tmp_path):
    converter = _absorbance_probe(tmp_path, params={'transmittance': True})
    assert converter.label['y'] == 'TRANSMITTANCE'
    # T = 10**(-A): the A = 2.0 band becomes 0.01, the 0.02 baseline ~0.955
    assert min(converter.ys) == pytest.approx(0.01, abs=1e-6)
    assert max(converter.ys) == pytest.approx(10 ** -0.02, abs=1e-6)


def test_both_instructions_compose(tmp_path):
    """Convert first, then mirror -- the order the client's words imply."""
    converter = _absorbance_probe(
        tmp_path, params={'transmittance': True, 'invert_y': True})
    assert converter.label['y'] == 'TRANSMITTANCE - inverted'


def test_transmittance_refuses_data_that_already_looks_like_transmittance(tmp_path):
    """The old heuristic, now a guard on an explicit request rather than a
    silent decision. A wrong call is a visible refusal, not corrupt output."""
    with pytest.raises(UnconvertibleSpectrum, match='already appears'):
        _probe(TRANSMITTANCE_SHAPED, 'INFRARED SPECTRUM', tmp_path,
               params={'transmittance': True})


def test_transmittance_refuses_data_that_is_not_absorbance_scaled(tmp_path):
    """JPK-948 declares ABSORBANCE but carries detector counts of 13-6168.
    10**(-6168) underflows, so converting would return a flat line labelled
    TRANSMITTANCE -- worse than anything this change fixes."""
    with pytest.raises(UnconvertibleSpectrum, match='not absorbance'):
        _probe(ABSORBANCE_SHAPED, 'INFRARED SPECTRUM', tmp_path,
               yunits='ABSORBANCE',
               params={'transmittance': True, 'jcamp_idx': 0})


# - - - provenance - - -

def test_flags_are_emitted_only_when_acted_on(tmp_path):
    plain = TechniqueComposer(_probe(
        ABSORBANCE_SHAPED, 'INFRARED SPECTRUM', tmp_path, yunits='ABSORBANCE'))
    assert '##$CSINVERTY' not in ''.join(plain.meta)
    assert '##$CSTRANSMITTANCE' not in ''.join(plain.meta)

    inverted = TechniqueComposer(_probe(
        ABSORBANCE_SHAPED, 'INFRARED SPECTRUM', tmp_path, yunits='ABSORBANCE',
        params={'invert_y': True}))
    assert '##$CSINVERTY=true' in ''.join(inverted.meta)


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
    assert _post(client, TRANSMITTANCE_SHAPED, invert_y='true').status_code == 200


@pytest.mark.parametrize('source, reason', [
    (TRANSMITTANCE_SHAPED, 'already appears'),
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
