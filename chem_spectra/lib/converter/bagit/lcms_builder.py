import logging
import os
import tempfile
from typing import List, Optional, Tuple

from chem_spectra.lib.composer.lcms_converter_app import (
    LCMSConverterAppComposer,
    has_uvvis_peak_marker,
)
from chem_spectra.lib.external.chemotion_converter_lcms import (
    lcms_preview_image_from_jdx_files,
)

logger = logging.getLogger(__name__)

_HEADER_SCAN_LINES = 200


def _read_jcamp_header(path: str) -> Tuple[Optional[str], Optional[str]]:
    data_type: Optional[str] = None
    scan_mode: Optional[str] = None
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as handle:
            for _ in range(_HEADER_SCAN_LINES):
                raw = handle.readline()
                if not raw:
                    break
                line = raw.strip()
                if not line.startswith('##'):
                    if line:
                        continue
                    break
                if '=' not in line:
                    continue
                key, value = line[2:].split('=', 1)
                key_norm = key.strip().lstrip('$').upper()
                value_norm = value.strip()
                if key_norm == 'DATA TYPE' and data_type is None:
                    data_type = value_norm
                elif key_norm == 'SCAN_MODE' and scan_mode is None:
                    scan_mode = value_norm
                if data_type and scan_mode:
                    break
    except OSError:
        return None, None
    return data_type, scan_mode


def _scan_mode_suffix(scan_mode: Optional[str]) -> Optional[str]:
    if not scan_mode:
        return None
    lowered = scan_mode.lower()
    if 'pos' in lowered:
        return 'pos'
    if 'neg' in lowered:
        return 'neg'
    return None


def _semantic_stem_for(data_type: Optional[str], scan_mode: Optional[str]) -> Optional[str]:
    if not data_type:
        return None
    dt = data_type.upper()
    if 'UV-VIS' in dt or 'UV/VIS' in dt or 'ULTRAVIOLET' in dt:
        return 'lcms_uvvis'
    polarity = _scan_mode_suffix(scan_mode)
    if 'MASS TIC' in dt or 'TOTAL ION CHROMATOGRAM' in dt or 'TOTAL ION CHROMATOGRAPHY' in dt:
        return f'lcms_tic_{polarity}' if polarity else 'lcms_tic'
    if 'MASS SPECTRUM' in dt:
        return f'lcms_mz_{polarity}' if polarity else 'lcms_mz'
    return None


def classify_lcms_stems(jcamp_paths: List[str]) -> List[str]:
    if not jcamp_paths:
        return []

    raw_stems = [os.path.splitext(os.path.basename(p))[0].replace('.', '_') for p in jcamp_paths]
    semantic_stems: List[Optional[str]] = []
    for path in jcamp_paths:
        if has_uvvis_peak_marker(path):
            semantic_stems.append('lcms_uvvis.peak')
            continue
        data_type, scan_mode = _read_jcamp_header(path)
        semantic_stems.append(_semantic_stem_for(data_type, scan_mode))

    counts: dict = {}
    for stem in semantic_stems:
        if stem:
            counts[stem] = counts.get(stem, 0) + 1

    resolved: List[str] = []
    for raw, semantic in zip(raw_stems, semantic_stems):
        if semantic and counts[semantic] == 1:
            resolved.append(semantic)
        else:
            resolved.append(raw)
    return resolved


def _copy_to_named_tmp(path: str) -> tempfile.NamedTemporaryFile:
    base = os.path.basename(path)
    suffix = '_{}'.format(base) if base else (os.path.splitext(path)[1] or '')
    tf = tempfile.NamedTemporaryFile(suffix=suffix)
    with open(path, 'rb') as src:
        tf.write(src.read())
    tf.seek(0)
    return tf


def _close_silently(handles):
    for handle in handles:
        try:
            handle.close()
        except Exception:
            pass


def build_lcms_composer(
    jcamp_paths,
    params=None,
) -> Optional[LCMSConverterAppComposer]:
    if not jcamp_paths:
        return None

    ordered_paths: List[str] = sorted(jcamp_paths)

    tmp_files: List[tempfile.NamedTemporaryFile] = []
    try:
        for path in ordered_paths:
            tmp_files.append(_copy_to_named_tmp(path))
    except Exception:
        _close_silently(tmp_files)
        return None

    try:
        preview = lcms_preview_image_from_jdx_files(tmp_files, params)
    except Exception:
        preview = None

    return LCMSConverterAppComposer(tmp_files, preview, params)


def append_lcms_group(
    lcms_paths,
    params,
    list_files,
    list_images,
    list_csv,
    list_composer,
    archive_stems=None,
):
    if not lcms_paths:
        return

    try:
        composer = build_lcms_composer(lcms_paths, params)
    except Exception as err:
        logger.warning(
            "Skip LC/MS group (%d files): builder raised %s",
            len(lcms_paths), err,
        )
        return

    if composer is None or not composer.data:
        logger.warning(
            "Skip LC/MS group (%d files): no usable composer",
            len(lcms_paths),
        )
        return

    try:
        composer.tf_jcamp()
    except Exception as err:
        logger.warning("LC/MS group: tf_jcamp failed (%s); using raw inputs", err)

    try:
        preview = composer.tf_img()
    except Exception:
        preview = None

    try:
        csv_file = composer.tf_csv()
    except Exception:
        csv_file = None

    files = list(composer.data or [])
    n = len(files)
    if n == 0:
        return

    stems = classify_lcms_stems([f.name for f in files])
    if archive_stems is not None:
        archive_stems.extend(stems)

    preview_idx = next(
        (i for i, s in enumerate(stems) if s == 'lcms_uvvis'),
        0,
    )
    image_slots = [None] * n
    image_slots[preview_idx] = preview

    list_files.extend(files)
    list_images.extend(image_slots)
    list_csv.extend([csv_file] + [None] * (n - 1))
    list_composer.extend([composer] * n)
