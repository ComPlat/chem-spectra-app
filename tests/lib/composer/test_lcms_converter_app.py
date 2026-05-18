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
