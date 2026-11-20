"""MS data tables are written as short as the values' own precision allows."""
import re

import numpy as np
import pytest

from chem_spectra.lib.shared.number_text import xy_lines

HALF_FLOAT32_STEP = 2.0 ** -24  # relative


def _parse(lines):
    return np.array([[float(v) for v in line.split(', ')] for line in lines])


def test_float32_data_is_written_short_and_reads_back_to_the_same_float32():
    rng = np.random.default_rng(0)
    xs = rng.uniform(50, 2000, 5000).astype(np.float32).astype(np.float64)
    ys = rng.uniform(0, 5e9, 5000).astype(np.float32).astype(np.float64)
    lines = xy_lines(xs, ys)
    back = _parse(lines)
    assert np.array_equal(back[:, 0].astype(np.float32), xs.astype(np.float32))
    assert np.array_equal(back[:, 1].astype(np.float32), ys.astype(np.float32))
    # as float64, within half a float32 step of the original
    assert np.max(np.abs(back[:, 0] - xs) / xs) <= HALF_FLOAT32_STEP
    assert sum(map(len, lines)) < sum(
        len('{}, {}\n'.format(x, y)) for x, y in zip(xs, ys)) * 0.75


def test_an_example():
    assert xy_lines([100.18072509765625], [81569.296875]) == [
        '100.180725, 81569.3\n']


def test_float64_data_is_written_as_before():
    xs = np.array([100.1807250976562, 0.1, 51.0120849609371])  # not float32
    ys = np.array([1.0, 2.5, 3.0])
    assert xy_lines(xs, ys) == [
        '{}, {}\n'.format(x, y) for x, y in zip(xs, ys)]


def test_one_float64_value_keeps_the_whole_column_as_before():
    xs = np.array([100.18072509765625, 0.1])  # 0.1 is not a float32
    assert xy_lines(xs, [1.0, 2.0])[0] == '100.18072509765625, 1.0\n'


@pytest.mark.parametrize('value, text', [
    (123456792.0, '123456790.0'),  # numpy would print 1.2345679e+08
    (1e-05, '1e-05'),              # repr's own notation, as before
    (2.0 ** 60, '1.1529215e+18'),  # beyond 1e16 repr uses e-notation too
])
def test_the_text_style_follows_repr(value, text):
    exact = float(np.float32(value))
    assert xy_lines([exact], [1.0])[0] == '{}, 1.0\n'.format(text)


def test_an_mzml_upload_shortens_only_what_is_float32(client):
    # svs813f1.mzML stores intensities as 32-bit floats and m/z as 64-bit, and
    # three m/z values per scan use that precision. The intensities come out
    # short; the m/z column keeps every digit.
    with open('./tests/fixtures/source/ms/svs813f1.mzML', 'rb') as f:
        response = client.post('/jcamp', data={'file': (f, 'svs813f1.mzML')},
                               content_type='multipart/form-data')
    assert response.status_code == 200
    points = [line.split(', ') for line in response.data.decode().splitlines()
              if re.match(r'^-?\d', line) and ', ' in line]
    assert points

    def digits(value):
        return len(re.sub(r'[^0-9]', '', value.split('e')[0]).strip('0'))

    assert max(digits(y) for _, y in points) <= 9
    assert '262.09028625488' in {x[:15] for x, _ in points}
