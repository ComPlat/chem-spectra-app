"""Tests for LC/MS handling inside BagIt archives.

These tests exercise the new ``build_lcms_composer`` integration in
``BagItBaseConverter``. They build self-contained BagIt directories on the fly
(rather than depending on binary fixtures) so that the layout under test is
explicit and the suite is robust to fixture changes.
"""
import os
import tempfile

import pytest

from chem_spectra.lib.composer.lcms_converter_app import LCMSConverterAppComposer
from chem_spectra.lib.converter.bagit.base import BagItBaseConverter
from chem_spectra.lib.converter.bagit.lcms_builder import build_lcms_composer


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
    """Materialize a minimal BagIt directory tree under ``root``.

    ``files`` maps filename → JCAMP content for ``data/``.
    """
    data_dir = os.path.join(root, 'data')
    os.makedirs(data_dir, exist_ok=True)
    for name, content in files.items():
        with open(os.path.join(data_dir, name), 'w', encoding='utf-8') as fh:
            fh.write(content)
    with open(os.path.join(root, 'bagit.txt'), 'w', encoding='utf-8') as fh:
        fh.write('BagIt-Version: 0.97\nTag-File-Character-Encoding: UTF-8\n')
    return root


def _write_flat_layout(root: str, files: dict) -> str:
    """Drop JCAMP files directly under ``root`` (no ``data/``, no ``bagit.txt``).

    Used to model a user upload that bundles ``.jdx`` files at the root of a
    plain zip archive (e.g. an LCMS export with ``uvvis.peak.jdx`` etc.).
    """
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
        # Image+csv attached only to the first item of the group.
        assert converter.images[1] is None and converter.images[2] is None
        assert converter.list_csv[0] is None
        assert converter.list_csv[1] is None and converter.list_csv[2] is None
        # Combined image must be None when LCMS is involved.
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
        assert converter.combined_image is None  # LCMS guard wins over CV combine.


def test_pure_cv_bagit_remains_strictly_master_compatible():
    """A BagIt with no LC/MS content must still go through CV combining and
    yield a combined image, exactly like on master."""
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
    # Each entry should be a NamedTemporaryFile-like handle pointing at a real
    # file on disk so that downstream consumers (zip writer, base64 encoder)
    # can re-read it.
    for handle in composer.data:
        assert hasattr(handle, 'name')
        assert os.path.isfile(handle.name)


def test_build_lcms_composer_skips_missing_path(tmp_path):
    missing = str(tmp_path / 'does_not_exist.jdx')
    composer = build_lcms_composer([missing], params=None)
    assert composer is None


def test_flat_layout_lcms_groups_into_single_composer():
    """A directory with LC/MS ``.jdx`` files at the root (no BagIt envelope)
    must be picked up by ``BagItBaseConverter`` and treated like a BagIt
    ``data/`` folder."""
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
        assert converter.combined_image is None  # LCMS guard


def test_flat_layout_ignores_non_jdx_siblings():
    """In flat layout we must not feed README/manifest files into the JCAMP
    parser; only ``.jdx`` files are picked up."""
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
    """UV/Vis traces use the LCMS composer path (with LC/MS/TIC), not NIComposer."""
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
    """Regression guard: presence of ``data/`` must keep the strict BagIt
    behaviour (root-level files are NOT picked up as fallback)."""
    with tempfile.TemporaryDirectory() as td:
        _write_bagit_layout(td, {
            'a_cv.jdx': CV_JDX_TEMPLATE.format(title='CV A'),
        })
        # Drop a stray .jdx at the root: it must be ignored because data/
        # exists and is the canonical source.
        with open(os.path.join(td, 'stray.jdx'), 'w', encoding='utf-8') as fh:
            fh.write(CV_JDX_TEMPLATE.format(title='Stray'))

        converter = BagItBaseConverter(td)

        assert converter.data is not None
        assert len(converter.data) == 1


def test_bagit_skips_lcms_group_when_builder_fails(monkeypatch):
    """If the LC/MS builder raises, the BagIt processing must continue with
    its remaining (non-LCMS) entries instead of crashing."""
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
        # Only the CV entry survives; LCMS group is skipped silently.
        assert converter.data is not None
        assert len(converter.data) == 1
