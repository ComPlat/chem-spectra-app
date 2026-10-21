"""A block whose data does not have the points it declares is refused.

JCAMP-DX uses the comma only as a delimiter. A file written with decimal
commas therefore parses without error into the wrong numbers -- every decimal
splits in two -- and was drawn and returned with a 200. What gives it away is
the count: across 235 conforming files the parsed count always equals
##NPOINTS, and a decimal comma breaks that every time.
"""

import io

import numpy as np
import pytest

from chem_spectra import create_app
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.data_parse import make_ms_data_xsys
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter
from chem_spectra.lib.converter.share import UnconvertibleSpectrum

IR_DECIMAL_HEADER = """##TITLE=decimal commas in the header
##JCAMP-DX=4.24
##DATA TYPE=INFRARED SPECTRUM
##XUNITS=1/CM
##YUNITS=ABSORBANCE
##XFACTOR=1,0
##YFACTOR=0,001
##FIRSTX=4000,0
##LASTX=3996,0
##DELTAX=-1,0
##NPOINTS=5
##XYDATA=(X++(Y..Y))
4000,0 120 130 140 150 160
##END=
"""

IR_DECIMAL_DATA = """##TITLE=decimal commas in the data
##JCAMP-DX=4.24
##DATA TYPE=INFRARED SPECTRUM
##XUNITS=1/CM
##YUNITS=ABSORBANCE
##XFACTOR=1
##YFACTOR=1
##FIRSTX=4000
##LASTX=3996
##DELTAX=-1
##NPOINTS=5
##XYDATA=(X++(Y..Y))
4000 0,120 0,130 0,140 0,150 0,160
##END=
"""

UV_SEMICOLON_PAIRS = """##TITLE=decimal commas, semicolon between x and y
##JCAMP-DX=5.00
##DATA TYPE=UV/VIS SPECTRUM
##XUNITS=NANOMETERS
##YUNITS=ABSORBANCE
##FIRSTX=400,0
##LASTX=404,0
##NPOINTS=5
##XYPOINTS=(XY..XY)
400,0;0,12
401,0;0,13
402,0;0,95
403,0;0,14
404,0;0,15
##END=
"""

MS_DECIMAL_PEAKS = """##TITLE=decimal commas in a peak table
##JCAMP-DX=5.00
##DATA TYPE=MASS SPECTRUM
##XUNITS=M/Z
##YUNITS=RELATIVE ABUNDANCE
##NPOINTS=2
##PEAK TABLE=(XY..XY)
15,5 420,25
16,0 300,5
##END=
"""

# The shape of jcamp-dx's doubleformat.jdx: decimal commas in two header LDRs,
# a conforming table underneath. The header pattern alone would refuse it.
MS_COMMA_HEADER_CONFORMING_DATA = """##TITLE=comma decimals in MAXX/MINX only
##JCAMP-DX=4.24
##DATA TYPE=MASS SPECTRUM
##XUNITS=M/Z
##YUNITS=RELATIVE ABUNDANCE
##FIRSTX=14
##LASTX=100
##MAXX=100,1
##MINX=14,1
##NPOINTS=4
##PEAK TABLE=(XY..XY)
14.1,4 15,31
16,1 100.1,1
##END=
"""

# What the 2018 fork's comma heuristic broke: integer pairs are pairs.
MS_INTEGER_PAIRS = """##TITLE=integer pairs
##JCAMP-DX=5.00
##DATA TYPE=MASS SPECTRUM
##XUNITS=M/Z
##YUNITS=RELATIVE ABUNDANCE
##NPOINTS=6
##PEAK TABLE=(XY..XY)
15,420 16,1201 26,58
27,513 28,1104 29,9999
##END=
"""


def _base(tmp_path, text, name='probe.jdx'):
    path = tmp_path / name
    path.write_text(text)
    return JcampBaseConverter(str(path), {})


@pytest.mark.parametrize('text, parsed, declared, comma', [
    (IR_DECIMAL_HEADER, 6, 5, '##FIRSTX=4000,0'),
    (IR_DECIMAL_DATA, 10, 5, None),
    (UV_SEMICOLON_PAIRS, 10, 5, '##FIRSTX=400,0'),
], ids=['ir-header', 'ir-data', 'uv-pairs'])
def test_a_spectrum_with_the_wrong_count_is_refused(tmp_path, text, parsed, declared, comma):
    with pytest.raises(UnconvertibleSpectrum) as refused:
        JcampTechniqueConverter(_base(tmp_path, text))
    reason = str(refused.value)
    assert '{} points'.format(parsed) in reason
    assert '##NPOINTS={}'.format(declared) in reason
    if comma:
        assert comma in reason and 'decimal commas' in reason
    else:
        assert 'decimal commas' not in reason


def test_a_mass_spectrum_with_the_wrong_count_is_refused(tmp_path):
    with pytest.raises(UnconvertibleSpectrum, match=r'4 points .* ##NPOINTS=2'):
        make_ms_data_xsys(_base(tmp_path, MS_DECIMAL_PEAKS))


def test_comma_decimals_in_the_header_alone_do_not_refuse(tmp_path):
    runs = make_ms_data_xsys(_base(tmp_path, MS_COMMA_HEADER_CONFORMING_DATA))
    np.testing.assert_allclose(runs[0], [[14.1, 4], [15, 31], [16, 1], [100.1, 1]])


def test_integer_pairs_stay_pairs(tmp_path):
    runs = make_ms_data_xsys(_base(tmp_path, MS_INTEGER_PAIRS))
    assert np.asarray(runs[0]).shape == (6, 2)
    np.testing.assert_allclose(runs[0][0], [15, 420])


def test_the_endpoint_says_why():
    client = create_app({'TESTING': True}).test_client()
    response = client.post(
        '/zip_jcamp_n_img',
        data={'file': (io.BytesIO(IR_DECIMAL_HEADER.encode()), 'ir.jdx')},
        content_type='multipart/form-data',
    )
    assert response.status_code == 422
    assert '##NPOINTS=5' in response.get_json()['error']
