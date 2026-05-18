import os
import tempfile

from chem_spectra.lib.composer.lcms_converter_app import LCMSConverterAppComposer
from chem_spectra.lib.converter.bagit.base import BagItBaseConverter
from chem_spectra.lib.converter.bagit.lcms_builder import (
    build_lcms_composer,
    classify_lcms_stems,
)


LCMS_JDX_TEMPLATE = """##TITLE={title}
##JCAMP-DX=5.00 $$ chemotion-converter-app (0.1.0)
##DATA TYPE=LC/MS
##DATA CLASS=XYPOINTS
##XUNITS=Time (s)
##YUNITS=Intensity
##NPOINTS=3
##XYPOINTS=(XY..XY)
1.0, 100.0
2.0, 200.0
3.0, 300.0
##END=
"""

CV_JDX_TEMPLATE = """##TITLE={title}
##JCAMP-DX=5.00 $$ chemotion-converter-app (0.1.0)
##DATA TYPE=CYCLIC VOLTAMMETRY
##DATA CLASS=XYPOINTS
##XUNITS=V vs Ref
##YUNITS=Ampere
##NPOINTS=3
##FIRSTX=0.0
##LASTX=2.0
##MINX=0.0
##MAXX=2.0
##MINY=0.0
##MAXY=2e-6
##FIRSTY=0.0
##XYPOINTS=(XY..XY)
0.0, 0.0
1.0, 1e-6
2.0, 2e-6
##END=
"""


def _write_bagit_layout(root: str, files: dict) -> str:
    data_dir = os.path.join(root, 'data')
    os.makedirs(data_dir, exist_ok=True)
    for name, content in files.items():
        with open(os.path.join(data_dir, name), 'w', encoding='utf-8') as fh:
            fh.write(content)
    with open(os.path.join(root, 'bagit.txt'), 'w', encoding='utf-8') as fh:
        fh.write('BagIt-Version: 0.97\nTag-File-Character-Encoding: UTF-8\n')
    return root


def _write_flat_layout(root: str, files: dict) -> str:
    for name, content in files.items():
        with open(os.path.join(root, name), 'w', encoding='utf-8') as fh:
            fh.write(content)
    return root


def test_lcms_only_bagit_groups_into_single_composer():
    with tempfile.TemporaryDirectory() as td:
        _write_bagit_layout(td, {
            'lc_a.jdx': LCMS_JDX_TEMPLATE.format(title='LCMS A'),
            'lc_b.jdx': LCMS_JDX_TEMPLATE.format(title='LCMS B'),
            'lc_c.jdx': LCMS_JDX_TEMPLATE.format(title='LCMS C'),
        })

        converter = BagItBaseConverter(td)

        assert converter.data is not None
        assert len(converter.data) == 3
        assert len(converter.images) == 3
        assert len(converter.list_csv) == 3
        assert converter.images[1] is None and converter.images[2] is None
        assert converter.list_csv[0] is None
        assert converter.list_csv[1] is None and converter.list_csv[2] is None
        assert converter.combined_image is None


def test_lcms_only_bagit_get_base64_data_returns_all_curves():
    with tempfile.TemporaryDirectory() as td:
        _write_bagit_layout(td, {
            'lc_a.jdx': LCMS_JDX_TEMPLATE.format(title='LCMS A'),
            'lc_b.jdx': LCMS_JDX_TEMPLATE.format(title='LCMS B'),
        })
        converter = BagItBaseConverter(td)
        base64_list = converter.get_base64_data()
        assert base64_list is not None
        assert len(base64_list) == 2


def test_mixed_bagit_keeps_existing_branches_and_adds_lcms_group():
    with tempfile.TemporaryDirectory() as td:
        _write_bagit_layout(td, {
            'a_cv.jdx': CV_JDX_TEMPLATE.format(title='CV A'),
            'b_cv.jdx': CV_JDX_TEMPLATE.format(title='CV B'),
            'c_lcms.jdx': LCMS_JDX_TEMPLATE.format(title='LCMS C'),
            'd_lcms.jdx': LCMS_JDX_TEMPLATE.format(title='LCMS D'),
        })

        converter = BagItBaseConverter(td)

        assert converter.data is not None
        assert len(converter.data) == 4
        assert converter.combined_image is None


def test_pure_cv_bagit_remains_strictly_master_compatible():
    with tempfile.TemporaryDirectory() as td:
        _write_bagit_layout(td, {
            'a_cv.jdx': CV_JDX_TEMPLATE.format(title='CV A'),
            'b_cv.jdx': CV_JDX_TEMPLATE.format(title='CV B'),
        })
        converter = BagItBaseConverter(td)
        assert len(converter.data) == 2
        assert converter.combined_image is not None


def test_build_lcms_composer_returns_none_on_empty_input():
    assert build_lcms_composer([], params=None) is None
    assert build_lcms_composer(None, params=None) is None


def test_build_lcms_composer_returns_composer_with_temp_handles(tmp_path):
    paths = []
    for name in ('lc_a.jdx', 'lc_b.jdx'):
        p = tmp_path / name
        p.write_text(LCMS_JDX_TEMPLATE.format(title=name))
        paths.append(str(p))

    composer = build_lcms_composer(paths, params=None)
    assert composer is not None
    assert isinstance(composer, LCMSConverterAppComposer)
    assert composer.data is not None
    assert len(composer.data) == 2
    for handle in composer.data:
        assert hasattr(handle, 'name')
        assert os.path.isfile(handle.name)


def test_build_lcms_composer_skips_missing_path(tmp_path):
    missing = str(tmp_path / 'does_not_exist.jdx')
    composer = build_lcms_composer([missing], params=None)
    assert composer is None


def test_flat_layout_lcms_groups_into_single_composer():
    with tempfile.TemporaryDirectory() as td:
        _write_flat_layout(td, {
            'tic.peak.jdx': LCMS_JDX_TEMPLATE.format(title='TIC peak'),
            'tic.edit.jdx': LCMS_JDX_TEMPLATE.format(title='TIC edit'),
            'ms.peak.jdx': LCMS_JDX_TEMPLATE.format(title='MS peak'),
            'ms.edit.jdx': LCMS_JDX_TEMPLATE.format(title='MS edit'),
        })

        converter = BagItBaseConverter(td)

        assert converter.data is not None
        assert len(converter.data) == 4
        assert converter.combined_image is None


def test_flat_layout_ignores_non_jdx_siblings():
    with tempfile.TemporaryDirectory() as td:
        _write_flat_layout(td, {
            'tic.peak.jdx': LCMS_JDX_TEMPLATE.format(title='TIC peak'),
            'README.txt': 'just some notes',
            'metadata.json': '{"foo": "bar"}',
        })

        converter = BagItBaseConverter(td)

        assert converter.data is not None
        assert len(converter.data) == 1


def test_flat_layout_uvvis_only_goes_through_lcms_group():
    uvvis_jdx = """##TITLE=UV-Vis only
##JCAMP-DX=5.00
##DATA TYPE=UV/VIS SPECTRUM
##DATA CLASS=XYPOINTS
##XUNITS=NANOMETERS
##YUNITS=ABSORBANCE
##NPOINTS=3
##XYPOINTS=(XY..XY)
200.0, 0.1
300.0, 0.2
400.0, 0.3
##END=
"""
    with tempfile.TemporaryDirectory() as td:
        _write_flat_layout(td, {'uvvis.peak.jdx': uvvis_jdx})

        converter = BagItBaseConverter(td)

        assert converter.data is not None
        assert len(converter.data) == 1
        assert converter.combined_image is None


def test_bagit_layout_unchanged_when_data_subdir_present():
    with tempfile.TemporaryDirectory() as td:
        _write_bagit_layout(td, {
            'a_cv.jdx': CV_JDX_TEMPLATE.format(title='CV A'),
        })
        with open(os.path.join(td, 'stray.jdx'), 'w', encoding='utf-8') as fh:
            fh.write(CV_JDX_TEMPLATE.format(title='Stray'))

        converter = BagItBaseConverter(td)

        assert converter.data is not None
        assert len(converter.data) == 1


def test_bagit_skips_lcms_group_when_builder_fails(monkeypatch):
    import chem_spectra.lib.converter.bagit.lcms_builder as lcms_builder_module

    def boom(_paths, _params):
        raise RuntimeError('synthetic LCMS failure')

    monkeypatch.setattr(lcms_builder_module, 'build_lcms_composer', boom)

    with tempfile.TemporaryDirectory() as td:
        _write_bagit_layout(td, {
            'a_cv.jdx': CV_JDX_TEMPLATE.format(title='CV A'),
            'b_lcms.jdx': LCMS_JDX_TEMPLATE.format(title='LCMS B'),
        })
        converter = BagItBaseConverter(td)
        assert converter.data is not None
        assert len(converter.data) == 1


UVVIS_NTUPLES = """##TITLE=Spectrum
##JCAMP-DX=5.00 $$ chemotion-converter-app (1.8.0)
##DATA TYPE=HPLC UV-VIS
##DATA CLASS=NTUPLES
##XUNITS=MINUTES
##YUNITS=SIGNAL
##NPOINTS=2
##PAGE=Wavelength= 210.0
##XYDATA=(XY..XY)
0.0, 1.0
1.0, 2.0
##END=
"""

MASS_TIC_POS = """##TITLE=Spectrum
##JCAMP-DX=5.00 $$ chemotion-converter-app (1.8.0)
##DATA TYPE=MASS TIC
##DATA CLASS=PEAK TABLE
##SCAN_MODE=positiv
##XUNITS=MINUTES
##YUNITS=COUNTS
##NPOINTS=2
##PEAK TABLE=(XY..XY)
0.0, 1.0
1.0, 2.0
##END=
"""

MASS_TIC_NEG = MASS_TIC_POS.replace('SCAN_MODE=positiv', 'SCAN_MODE=negativ')

MASS_SPEC_POS = """##TITLE=Spectrum
##JCAMP-DX=5.00 $$ chemotion-converter-app (1.8.0)
##DATA TYPE=MASS SPECTRUM
##DATA CLASS=NTUPLES
##SCAN_MODE=positiv
##XUNITS=m/z
##YUNITS=Intensity
##NPOINTS=2
##PEAK TABLE=(XY..XY)
100.0, 10.0
200.0, 20.0
##END=
"""

MASS_SPEC_NEG = MASS_SPEC_POS.replace('SCAN_MODE=positiv', 'SCAN_MODE=negativ')

UNKNOWN_TYPE = """##TITLE=Spectrum
##JCAMP-DX=5.00
##DATA TYPE=SOMETHING NEW
##NPOINTS=1
##XYDATA=(XY..XY)
0.0, 0.0
##END=
"""


def _write_jdx(directory, name, content):
    path = os.path.join(directory, name)
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(content)
    return path


def test_classify_lcms_stems_returns_empty_list_for_no_inputs():
    assert classify_lcms_stems([]) == []


def test_classify_lcms_stems_maps_openlab_layout_to_eln_recognisable_stems(tmp_path):
    paths = [
        _write_jdx(tmp_path, '01.jdx', MASS_TIC_NEG),
        _write_jdx(tmp_path, '02.jdx', MASS_TIC_POS),
        _write_jdx(tmp_path, 'NTUPLES0.jdx', UVVIS_NTUPLES),
        _write_jdx(tmp_path, 'NTUPLES1.jdx', MASS_SPEC_NEG),
        _write_jdx(tmp_path, 'NTUPLES2.jdx', MASS_SPEC_POS),
    ]
    assert classify_lcms_stems(paths) == [
        'lcms_tic_neg',
        'lcms_tic_pos',
        'lcms_uvvis',
        'lcms_mz_neg',
        'lcms_mz_pos',
    ]


def test_classify_lcms_stems_keeps_source_name_for_unknown_type(tmp_path):
    path = _write_jdx(tmp_path, 'mystery.jdx', UNKNOWN_TYPE)
    assert classify_lcms_stems([path]) == ['mystery']


def test_classify_lcms_stems_keeps_source_name_when_scan_mode_missing(tmp_path):
    content = MASS_TIC_POS.replace('##SCAN_MODE=positiv\n', '')
    path = _write_jdx(tmp_path, 'tic.jdx', content)
    assert classify_lcms_stems([path]) == ['lcms_tic']


def test_classify_lcms_stems_keeps_source_name_on_semantic_collision(tmp_path):
    paths = [
        _write_jdx(tmp_path, 'first.jdx', MASS_TIC_POS),
        _write_jdx(tmp_path, 'second.jdx', MASS_TIC_POS),
    ]
    assert classify_lcms_stems(paths) == ['first', 'second']


def test_classify_lcms_stems_handles_unreadable_path_gracefully(tmp_path):
    paths = [
        str(tmp_path / 'does_not_exist.jdx'),
        _write_jdx(tmp_path, 'uvvis.jdx', UVVIS_NTUPLES),
    ]
    assert classify_lcms_stems(paths) == ['does_not_exist', 'lcms_uvvis']


def test_bagit_archive_stems_include_both_uvvis_source_and_prebaked_peak(tmp_path):
    data_dir = tmp_path / 'data'
    data_dir.mkdir()

    _write_jdx(data_dir, '01.jdx', MASS_TIC_NEG)
    _write_jdx(data_dir, 'NTUPLES0.jdx', UVVIS_NTUPLES)
    _write_jdx(data_dir, 'NTUPLES1.jdx', MASS_SPEC_POS)
    with open(tmp_path / 'bagit.txt', 'w', encoding='utf-8') as fh:
        fh.write('BagIt-Version: 0.97\n')

    converter = BagItBaseConverter(str(tmp_path))
    assert converter.data is not None
    assert converter.archive_entry_stems == [
        'lcms_tic_neg',
        'lcms_uvvis',
        'lcms_uvvis.peak',
        'lcms_mz_pos',
    ]
    assert converter.combined_image is None
    assert converter.images[0] is None
    assert converter.images[1] is not None
    assert converter.images[2] is None
    assert converter.images[3] is None


def test_classify_lcms_stems_detects_re_uploaded_peak_via_marker(tmp_path):
    peak_content = """##TITLE=Spectrum
##JCAMP-DX=5.00
##DATA TYPE=LC/MS
##DATA CLASS=PEAK TABLE

$$ === CHEMSPECTRA UVVIS PEAK TABLE ===
##PAGE=210.0
##DATA TABLE= (XY..XY), PEAKS
0.0, 1.0;
1.0, 2.0;
##END=
"""
    path = _write_jdx(tmp_path, 'cf9d2f40_table_lcms_uvvis_peak.jdx', peak_content)
    assert classify_lcms_stems([path]) == ['lcms_uvvis.peak']
