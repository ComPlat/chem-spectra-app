"""Integral and multiplet records written by MestReNova, not by ChemSpectra.

Mnova's `$OBSERVED MULTIPLETS` and `$OBSERVED MULTIPLETS PEAKS` have the same
columns ChemSpectra writes, but each opens with a column-header line --
`(I, X, Y, Hs, Vs, An, H, C, N, (J))` and `(I, X, Y)` -- which ChemSpectra's
own records do not. The header was read as data: `float('X')`, an HTML 500
on every upload of such a file. The fixture is synthetic: MestReNova's record
layout (header lines, `$$` column notes, an AREAS record, a LINK with a peak
table block) around made-up data.
"""
import base64
import io
import re

from chem_spectra.lib.converter.jcamp.technique import peak_table_rows

FIXTURE = './tests/fixtures/source/mnova/mnova_1h_multiplets.jdx'


def _post(client, url):
    with open(FIXTURE, 'rb') as f:
        return client.post(url, data={'file': (f, 'mnova_1h_multiplets.jdx')},
                           content_type='multipart/form-data')


def _record(jcamp, label):
    """The value lines of `##<label>=`, up to the next record or comment."""
    lines = jcamp.split('\n')
    start = next(i for i, line in enumerate(lines)
                 if line.startswith('##{}='.format(label)))
    body = []
    for line in lines[start + 1:]:
        if line.startswith(('##', '$$')):
            break
        body.append(line)
    return [line for line in body if line.strip()]


def test_file_convert_answers_with_the_spectrum(client):
    response = _post(client, '/api/v1/chemspectra/file/convert')
    assert response.status_code == 200
    jcamp = base64.b64decode(response.json['jcamp']).decode()
    assert '##DATA TYPE=NMR SPECTRUM' in jcamp


def test_the_multiplets_come_back_as_data_rows_only(client):
    response = _post(client, '/jcamp')
    assert response.status_code == 200
    jcamp = response.data.decode()
    multiplets = _record(jcamp, '$OBSERVEDMULTIPLETS')
    peaks = _record(jcamp, '$OBSERVEDMULTIPLETSPEAKS')
    assert len(multiplets) == 3
    assert multiplets[0].startswith('(1, 7.46')
    assert all(re.match(r'^\(\d', row) for row in multiplets + peaks)
    # every multiplet keeps its peaks, under its own index
    assert {row.split(',')[0] for row in peaks} == {'(1', '(2', '(3'}


def test_the_image_draws_the_multiplets(client):
    response = _post(client, '/zip_jcamp_n_img')
    assert response.status_code == 200


def test_rows_are_recognised_by_shape():
    record = '(I, X, Y)\n(1, 7.90, 326.9)\n\n(2, 7.83, 315.9)\n'
    assert peak_table_rows(record) == '(1, 7.90, 326.9)\n(2, 7.83, 315.9)'


def test_an_integral_header_is_dropped_and_a_missing_one_costs_no_row():
    with_header = '(X Y Z)\n(0.23, 0.31, 18)\n(3.90, 3.98, 3.02)'
    without = '(0.23, 0.31, 18)\n(3.90, 3.98, 3.02)'
    assert peak_table_rows(with_header) == peak_table_rows(without) == (
        '(0.23, 0.31, 18)\n(3.90, 3.98, 3.02)')


def test_a_name_with_an_escaped_comma_stays_one_field(client):
    jcamp = _post(client, '/jcamp').data.decode()
    first = _record(jcamp, '$OBSERVEDMULTIPLETS')[0]
    assert first.endswith("t, 2'&comma;6', 7.5000)")
