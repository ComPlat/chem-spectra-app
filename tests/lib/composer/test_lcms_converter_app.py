import os
import tempfile

import pytest

from chem_spectra.lib.composer.lcms_converter_app import (
    LCMSConverterAppComposer,
    has_uvvis_peak_marker,
)


UVVIS_NTUPLES_JDX = """##TITLE=Spectrum
##JCAMP-DX=5.00 $$ chemotion-converter-app (1.8.0)
##DATA TYPE=HPLC UV-VIS
##DATA CLASS=NTUPLES
##ORIGIN=
##OWNER=
##XUNITS=MINUTES
##YUNITS=SIGNAL
##ALLWAVES=[210.0, 254.0]
##NTUPLES=MULTIDIMENSIONAL
##VAR_NAME=Wavelength, MINUTES, SIGNAL
##SYMBOL=T, X, Y
##VAR_TYPE=PAGE, X, Y
##VAR_DIM=2, 3, 3
##UNITS=, MINUTES, SIGNAL
##PAGE=Wavelength= 210.0
##XYDATA=(XY..XY)
0.0, 10.0
1.0, 12.0
2.0, 11.0
##PAGE=Wavelength= 254.0
##XYDATA=(XY..XY)
0.0, 5.0
1.0, 6.5
2.0, 6.0
##END=
"""

MASS_TIC_JDX = """##TITLE=Spectrum
##JCAMP-DX=5.00 $$ chemotion-converter-app (1.8.0)
##DATA TYPE=MASS TIC
##DATA CLASS=PEAK TABLE
##XUNITS=MINUTES
##YUNITS=COUNTS
##SCAN_MODE=positiv
##NPOINTS=3
##PEAK TABLE=(XY..XY)
0.0, 100.0
1.0, 200.0
2.0, 150.0
##END=
"""

# jcampconverter (react-spectra-editor's frontend JCAMP parser) can't parse
# the ##XYPOINTS= LDR — only ##XYDATA=/##DATA TABLE=/##PEAK TABLE=. Some
# chemotion-converter-app TIC exports use ##XYPOINTS= as the primary data
# block, which this composer must normalize since it otherwise passes
# source bytes through untouched.
XYPOINTS_TIC_JDX = """##TITLE=Spectrum
##JCAMP-DX=5.00 $$ chemotion-converter-app (1.9.2)
##DATA TYPE=MASS TIC
##DATA CLASS=XYPOINTS
##XUNITS=MINUTES
##YUNITS=ARBITRARY UNITS
##FIRSTX=0.0
##LASTX=2.0
##NPOINTS=3
##XYPOINTS= (XY..XY)
0.0, 100;
1.0, 5000;
2.0, 200;
##END=
"""


def _named_tmp(content: str) -> tempfile.NamedTemporaryFile:
    tf = tempfile.NamedTemporaryFile(suffix=".jdx", delete=False)
    tf.write(content.encode("utf-8"))
    tf.seek(0)
    return tf


def _file_contains_marker(jdx_file) -> bool:
    return has_uvvis_peak_marker(jdx_file.name)


@pytest.fixture
def uvvis_tmp():
    tf = _named_tmp(UVVIS_NTUPLES_JDX)
    yield tf
    try:
        tf.close()
        if os.path.exists(tf.name):
            os.unlink(tf.name)
    except OSError:
        pass


@pytest.fixture
def mass_tic_tmp():
    tf = _named_tmp(MASS_TIC_JDX)
    yield tf
    try:
        tf.close()
        if os.path.exists(tf.name):
            os.unlink(tf.name)
    except OSError:
        pass


@pytest.fixture
def xypoints_tic_tmp():
    tf = _named_tmp(XYPOINTS_TIC_JDX)
    yield tf
    try:
        tf.close()
        if os.path.exists(tf.name):
            os.unlink(tf.name)
    except OSError:
        pass


def test_init_prebakes_uvvis_peak_alongside_ntuples_source(uvvis_tmp):
    composer = LCMSConverterAppComposer([uvvis_tmp], None, None)

    assert len(composer.data) == 2
    source, peak = composer.data
    assert source is uvvis_tmp
    assert not _file_contains_marker(source)
    assert peak is not uvvis_tmp
    assert peak.name.lower().endswith("peak.jdx")
    assert _file_contains_marker(peak)


def test_init_is_idempotent_when_peak_marker_already_present(uvvis_tmp):
    composer = LCMSConverterAppComposer([uvvis_tmp], None, None)
    peak = composer.data[1]
    peak_path = peak.name
    peak_size = os.path.getsize(peak_path)

    with open(peak_path, "rb") as src:
        payload = src.read()
    fresh = tempfile.NamedTemporaryFile(suffix=".jdx", delete=False)
    fresh.write(payload)
    fresh.seek(0)
    try:
        second = LCMSConverterAppComposer([fresh], None, None)
        assert len(second.data) == 1
        assert second.data[0] is fresh
        assert os.path.getsize(second.data[0].name) == peak_size
    finally:
        fresh.close()
        if os.path.exists(fresh.name):
            os.unlink(fresh.name)


def test_init_does_not_misdetect_mass_tic_as_uvvis(mass_tic_tmp):
    composer = LCMSConverterAppComposer([mass_tic_tmp], None, None)

    assert len(composer.data) == 1
    kept = composer.data[0]
    assert kept is mass_tic_tmp
    assert not _file_contains_marker(kept)


def test_init_handles_mixed_batch_inserts_peak_next_to_uvvis(uvvis_tmp, mass_tic_tmp):
    composer = LCMSConverterAppComposer([mass_tic_tmp, uvvis_tmp], None, None)

    assert len(composer.data) == 3
    assert composer.data[0] is mass_tic_tmp
    assert composer.data[1] is uvvis_tmp
    assert composer.data[2] is not uvvis_tmp
    assert composer.data[2].name.lower().endswith("peak.jdx")
    assert _file_contains_marker(composer.data[2])


def test_init_prebakes_uvvis_despite_non_uvvis_peak_sibling(uvvis_tmp, mass_tic_tmp):
    tic_peak = _named_tmp(MASS_TIC_JDX)
    tic_peak.close()
    os.rename(tic_peak.name, tic_peak.name.replace(".jdx", ".peak.jdx"))
    tic_peak.name = tic_peak.name.replace(".jdx", ".peak.jdx")
    try:
        composer = LCMSConverterAppComposer([tic_peak, uvvis_tmp], None, None)

        assert len(composer.data) == 3
        assert composer.data[0] is tic_peak
        assert composer.data[1] is uvvis_tmp
        assert composer.data[2].name.lower().endswith("peak.jdx")
        assert _file_contains_marker(composer.data[2])
    finally:
        if os.path.exists(tic_peak.name):
            os.unlink(tic_peak.name)


def test_preview_pipeline_can_now_extract_uvvis(uvvis_tmp):
    from chem_spectra.lib.external.chemotion_converter_lcms import (
        _extract_uvvis_from_peak_content,
    )

    composer = LCMSConverterAppComposer([uvvis_tmp], None, None)
    with open(composer.data[1].name, "r", encoding="utf-8", errors="ignore") as h:
        content = h.read()

    result = _extract_uvvis_from_peak_content(content)
    assert result is not None
    xs, ys, _edit_peaks, _integrations, _wl = result
    assert len(xs) > 0 and len(ys) > 0


def test_init_normalizes_xypoints_ldr_to_xydata(xypoints_tic_tmp):
    # jcampconverter (the frontend parser) can't read ##XYPOINTS=. This
    # composer passes source bytes through untouched otherwise, so it must
    # rewrite the LDR (and the DATA CLASS declaring it) to ##XYDATA=,
    # matching what NIComposer/MSComposer already normalize to.
    LCMSConverterAppComposer([xypoints_tic_tmp], None, None)

    with open(xypoints_tic_tmp.name, "r", encoding="utf-8", errors="ignore") as h:
        content = h.read()

    assert "##XYPOINTS=" not in content
    assert "##DATA CLASS=XYDATA" in content
    assert "##XYDATA= (XY..XY)" in content
    # data values themselves must be untouched
    assert "0.0, 100;" in content
    assert "1.0, 5000;" in content
    assert "2.0, 200;" in content


def test_init_leaves_non_xypoints_ldr_untouched(mass_tic_tmp):
    with open(mass_tic_tmp.name, "r", encoding="utf-8", errors="ignore") as h:
        before = h.read()

    LCMSConverterAppComposer([mass_tic_tmp], None, None)

    with open(mass_tic_tmp.name, "r", encoding="utf-8", errors="ignore") as h:
        after = h.read()

    assert before == after


def test_init_does_not_touch_triple_hash_xypoints_metadata_echo():
    # ###XYPOINTS= (three hashes) is chem-spectra-app's own original-metadata
    # echo format, not a live LDR jcampconverter would try to parse — the
    # normalization must only match the real double-hash ##XYPOINTS= LDR.
    content = (
        "##TITLE=Spectrum\n"
        "##DATA TYPE=CYCLIC VOLTAMMETRY\n"
        "##DATA CLASS=XYDATA\n"
        "##XYDATA= (X++(Y..Y))\n"
        "1 2 3\n"
        "$$ === CHEMSPECTRA ORIGINAL METADATA ===\n"
        "###XYPOINTS= (XY..XY)\n"
        "0.0, 100;\n"
        "##END=\n"
    )
    tf = _named_tmp(content)
    try:
        LCMSConverterAppComposer([tf], None, None)
        with open(tf.name, "r", encoding="utf-8", errors="ignore") as h:
            after = h.read()
        assert after == content
    finally:
        tf.close()
        if os.path.exists(tf.name):
            os.unlink(tf.name)
