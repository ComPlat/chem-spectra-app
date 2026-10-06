"""A JCAMP holding more than one dimension is refused, not flattened.

ChemSpectra reads one row of a 2D file and returns it as a 1D spectrum, with
the acquisition time axis labelled as the spectrum's. Measured on
`origin/master` before this change: a synthetic `nD NMR FID` with
`##NUM DIM= 2` returned **200** from `/api/v1/chemspectra/file/convert`,
`/zip_jcamp` and `/zip_jcamp_n_img`, the composed file carrying
`##DATA TYPE=ND NMR FID` and `##XUNITS=SECONDS`.

Nothing in the answer says the result is meaningless, which is the whole
problem: a caller cannot tell a flattened 2D file from a 1D one.

The header is read as text, before nmrglue parses it. Not because nmrglue
loses the records -- it keeps `NUMDIM` and `DATATYPE` on a small file -- but
because a real 2D dataset is large and need not parse at all, and there is no
reason to spend that read on a file that is going to be refused.
"""

import io

import pytest

from chem_spectra.lib.converter.jcamp.base import (
    JcampBaseConverter, declared_dimensions, declares_nmr,
)
from chem_spectra.lib.converter.jcamp.technique import UnconvertibleSpectrum

ONE_D = './tests/fixtures/source/1H.dx'                 # no ##NUM DIM at all
ONE_D_DECLARED = './tests/fixtures/source/mnova/STM212_H.jcamp'   # NUMDIM<tab>1


def _two_d(tmp_path, name='twod.dx', datatype='nD NMR FID',
           num_dim='##NUM DIM= 2'):
    ys = ' '.join(str(v) for v in [100, 200, 400, 300, 150, 120, 110, 105])
    lines = ['##TITLE=synthetic 2D', '##JCAMP-DX=6.00',
             '##DATA TYPE= ' + datatype]
    if num_dim:
        lines.append(num_dim)
    lines += ['##DATA CLASS= NTUPLES', '##XUNITS=SECONDS',
              '##YUNITS=ARBITRARY', '##FIRSTX=0', '##LASTX=7', '##NPOINTS=8',
              '##FIRSTY=100', '##XYDATA=(X++(Y..Y))', '0 ' + ys, '##END=']
    target = tmp_path / name
    target.write_text('\n'.join(lines) + '\n')
    return target


# - - - the detection function, on its own - - -

@pytest.mark.parametrize('header, expected', [
    ('##NUM DIM= 2', 2),
    ('##NUMDIM=\t2', 2),
    ('##NUM_DIM = 3', 3),
    ('##num dim=2', 2),
    ('##NUM DIM= 1', 1),
    ('##NUMDIM=\t1', 1),
    ('', None),
])
def test_the_declared_dimension_count_is_read_from_the_label(header, expected):
    """JCAMP-DX ignores spaces, dashes, underscores and slashes inside a
    label, so `NUM DIM`, `NUMDIM` and `NUM_DIM` are one record. Five fixtures
    here write it `##NUMDIM=<tab>1`."""
    text = '##TITLE=t\n{}\n##END=\n'.format(header) if header else '##TITLE=t\n'
    assert declared_dimensions(text) == expected


@pytest.mark.parametrize('datatype', [
    'nD NMR FID', 'nD NMR SPECTRUM', '2D NMR SPECTRUM', 'ND NMR FID'])
def test_an_nd_datatype_counts_as_multi_dimensional(datatype):
    """`nD` is the JCAMP-DX 6 spelling; some vendors write the number."""
    text = '##TITLE=t\n##DATA TYPE= {}\n##END=\n'.format(datatype)
    assert (declared_dimensions(text) or 0) > 1


@pytest.mark.parametrize('datatype', ['NMR SPECTRUM', 'NMR FID', 'INFRARED SPECTRUM'])
def test_a_one_dimensional_datatype_does_not(datatype):
    text = '##TITLE=t\n##DATA TYPE= {}\n##END=\n'.format(datatype)
    assert declared_dimensions(text) in (None, 1)


# - - - the converter refuses - - -

@pytest.mark.parametrize('kwargs', [
    {'datatype': 'nD NMR FID'},
    {'datatype': 'nD NMR SPECTRUM'},
    {'datatype': 'NMR SPECTRUM', 'num_dim': '##NUMDIM=\t2'},
    {'datatype': '2D NMR SPECTRUM', 'num_dim': None},
])
def test_a_multi_dimensional_file_is_refused(kwargs, tmp_path):
    source = _two_d(tmp_path, **kwargs)
    with pytest.raises(UnconvertibleSpectrum, match='2D'):
        JcampBaseConverter(str(source))


@pytest.mark.parametrize('source', [ONE_D, ONE_D_DECLARED])
def test_one_dimensional_files_are_untouched(source):
    """`STM212_H.jcamp` declares `##NUMDIM=<tab>1`, so the rule has to parse
    the number rather than react to the record being present."""
    assert JcampBaseConverter(source).typ


# - - - and every endpoint says so - - -

def _post(client, route, path, name=None):
    import os
    with open(path, 'rb') as handle:
        data = {'file': (io.BytesIO(handle.read()),
                         name or os.path.basename(path))}
    return client.post(route, content_type='multipart/form-data', data=data)


@pytest.mark.parametrize('route', [
    '/api/v1/chemspectra/file/convert',
    '/zip_jcamp',
    '/zip_jcamp_n_img',
])
def test_every_endpoint_refuses_with_a_reason(route, client, tmp_path):
    """422 with the reason in `error`, which is what the ELN and the
    standalone client show to the user. All three answered 200 before."""
    import json
    response = _post(client, route, str(_two_d(tmp_path)))
    assert response.status_code == 422
    assert '2D' in json.loads(response.data)['error']


@pytest.mark.parametrize('route', [
    '/api/v1/chemspectra/file/save',
    '/api/v1/chemspectra/file/refresh',
])
def test_the_save_paths_refuse_it_too(route, client, tmp_path):
    """Saving and refreshing re-read the uploaded file, so they reach the
    same guard. The standalone client downloads a save response as a zip, so
    a 200 carrying a flattened spectrum would be written to disk.

    These two take the file as `dst_list`, with `src` alongside on save.
    """
    import json
    body = _two_d(tmp_path).read_bytes()
    data = {'dst_list': (io.BytesIO(body), 'twod.dx')}
    if route.endswith('/save'):
        data['src'] = (io.BytesIO(body), 'twod.dx')
    response = client.post(route, content_type='multipart/form-data',
                           data=data)
    assert response.status_code == 422
    assert '2D' in json.loads(response.data)['error']


def test_an_archive_carrying_one_is_refused_whole(client, tmp_path):
    """One 2D member refuses the upload, rather than a partial result that
    looks complete. The same call this repository already makes for an upload
    holding two BagIt archives: silently processing part of it is the defect,
    not the remedy.
    """
    import json
    import zipfile
    archive = tmp_path / 'mixed.zip'
    with zipfile.ZipFile(archive, 'w') as zf:
        zf.write('./tests/fixtures/source/1H.dx', 'data/one.jdx')
        zf.writestr('data/two.jdx', _two_d(tmp_path).read_text())
        zf.writestr('bagit.txt', 'BagIt-Version: 0.97\n')
    response = _post(client, '/zip_jcamp_n_img', str(archive))
    assert response.status_code == 422
    assert '2D' in json.loads(response.data)['error']


# - - - agreeing with the ELN's rule, and where we deliberately differ - - -

@pytest.mark.parametrize('ending', ['\n', '\r\n', '\r'])
def test_every_line_ending_is_read(ending, tmp_path):
    """A CR-only file -- classic Mac, and some instrument exports -- is one
    long line to a multiline regex, so every record after the first was
    invisible and the file was accepted."""
    target = tmp_path / 'endings.dx'
    target.write_bytes(ending.join([
        '##TITLE=t', '##DATA TYPE= NMR SPECTRUM', '##NUM DIM= 2',
        '##XUNITS=SECONDS', '##YUNITS=ARBITRARY', '##NPOINTS=4',
        '##XYDATA=(X++(Y..Y))', '0 1 2 3 4', '##END=', '']).encode())
    with pytest.raises(UnconvertibleSpectrum, match='2D'):
        JcampBaseConverter(str(target))


def test_a_multi_dimensional_non_nmr_file_is_refused_without_naming_nmr(
        tmp_path):
    """chemotion_ELN requires NMR before it diverts a file, because it is
    choosing between ChemSpectra and NMRium. The question here is different
    -- whether this app can represent the data at all -- and it cannot, for
    any technique: `xs`/`ys` hold one curve.

    So a multi-dimensional UV/VIS file is still refused, but the reason does
    not call it NMR and does not send the user to NMRium, which would not
    read it either.
    """
    source = _two_d(tmp_path, datatype='UV/VIS SPECTRUM',
                    num_dim='##NUM DIM= 2')
    with pytest.raises(UnconvertibleSpectrum, match='declares 2 dimensions'):
        JcampBaseConverter(str(source))


def test_an_observed_nucleus_is_enough_to_name_it_nmr(tmp_path):
    """The datatype need not say NMR: `##.OBSERVE NUCLEUS=` does."""
    target = tmp_path / 'nucleus.dx'
    target.write_text('\n'.join([
        '##TITLE=t', '##DATA TYPE= SPECTRUM', '##.OBSERVE NUCLEUS= ^1H',
        '##NUM DIM= 2', '##XUNITS=SECONDS', '##YUNITS=ARBITRARY',
        '##NPOINTS=4', '##XYDATA=(X++(Y..Y))', '0 1 2 3 4', '##END=', '']))
    with pytest.raises(UnconvertibleSpectrum, match='NMRium'):
        JcampBaseConverter(str(target))


# - - - the label rule, as JCAMP-DX states it - - -

@pytest.mark.parametrize('label', [
    '##NUM DIM', '##NUMDIM', '##NUM_DIM', '##NUM-DIM',
    '##N-UM D/IM', '##N U M D I M', '##num dim',
])
def test_a_label_is_compared_with_its_separators_removed(label):
    """4.24 (5.1) compares labels with spaces, dashes, underscores and
    slashes removed, and without regard to case -- anywhere in the label, not
    only between its words. Matching each spelling with an optional separator
    between `NUM` and `DIM` left `##N-UM D/IM=` through, which is the same
    record.
    """
    assert declared_dimensions('##TITLE=t\n{}= 2\n'.format(label)) == 2


def test_the_datatype_says_how_many_when_it_knows():
    """`nD` does not say which n, so two is the least it can mean; `3D` does
    say."""
    def dims(datatype):
        return declared_dimensions(
            '##TITLE=t\n##DATA TYPE= {}\n'.format(datatype))
    assert dims('nD NMR FID') == 2
    assert dims('2D NMR SPECTRUM') == 2
    assert dims('3D NMR FID') == 3
    assert dims('NMR SPECTRUM') is None


@pytest.mark.parametrize('header, expected', [
    ('##DATA TYPE= NMR SPECTRUM', True),
    ('##DATA TYPE= nD NMR FID', True),
    ('##.OBSERVE NUCLEUS= ^1H', True),
    ('##.OBSERVE-NUCLEUS= ^1H', True),
    ('##DATA TYPE= UV/VIS SPECTRUM', False),
    ('##DATA TYPE= MASS SPECTRUM', False),
])
def test_what_counts_as_nmr_for_the_wording(header, expected):
    assert declares_nmr('##TITLE=t\n{}\n'.format(header)) is expected


# - - - a refusal must not leave anything behind on the shared figure - - -

def test_a_refused_overlay_leaves_the_figure_clean(client, tmp_path):
    """`/combine_images` plots onto the module-global pyplot figure.

    The converter refuses a 2D member from inside that loop, after the
    earlier files have already been drawn, and nothing on the raising path
    clears the figure. The curves stayed on `plt.gca()` and were drawn into
    the *next* image the worker rendered -- a different request, whose
    attachment image the ELN then stores. One request corrupting another's
    saved preview is worse than the flattening this PR set out to fix.

    The files are named so the 1D one sorts first and is plotted before the
    2D one raises.
    """
    import json

    import matplotlib.pyplot as plt

    with open(ONE_D, 'rb') as handle:
        one_d = handle.read()
    two_d = _two_d(tmp_path).read_bytes()

    plt.clf()
    plt.cla()
    response = client.post(
        '/combine_images',
        content_type='multipart/form-data',
        data={'files[]': [(io.BytesIO(one_d), 'a_1d.dx'),
                          (io.BytesIO(two_d), 'b_2d.dx')]},
    )
    assert response.status_code == 422
    assert '2D' in json.loads(response.data)['error']
    # the refusal names the member, because the caller sent several files
    assert 'b_2d' in json.loads(response.data)['error']
    assert len(plt.gca().lines) == 0, (
        'a refused overlay left %d curve(s) on the shared figure'
        % len(plt.gca().lines))


def test_a_refused_overlay_does_not_convert_the_other_members(client,
                                                              tmp_path):
    """And it refuses before any of them is read.

    Every member used to be parsed by nmrglue, converted and rendered to a
    3200x1800 PNG before the 2D one raised: a twenty-member upload paid
    nineteen full conversions for a 422. The header check costs a read of the
    first 64 KB.
    """
    import json

    unparsable = io.BytesIO(b'not a jcamp at all\n')
    two_d = _two_d(tmp_path).read_bytes()
    response = client.post(
        '/combine_images',
        content_type='multipart/form-data',
        data={'files[]': [(unparsable, 'a_broken.dx'),
                          (io.BytesIO(two_d), 'b_2d.dx')]},
    )
    # the unreadable member never gets as far as being read: the answer is
    # about the dimensions, not about the broken file
    assert response.status_code == 422
    assert '2D' in json.loads(response.data)['error']
