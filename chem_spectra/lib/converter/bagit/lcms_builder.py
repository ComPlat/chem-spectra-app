"""LCMS support for BagIt datasets.

Keeps the multi-file LC/MS handling out of ``BagItBaseConverter`` so the
generic NI/MS pipeline stays untouched. Two public entry points:

- ``build_lcms_composer``: turn a list of LC/MS JCAMP paths into a single
  :class:`LCMSConverterAppComposer` (with its preview).
- ``append_lcms_group``: glue helper used by ``BagItBaseConverter.__read``
  to extend the per-file ``list_files / list_images / list_csv /
  list_composer`` lists with a freshly-built composer, with safe fallbacks.
"""

import logging
import os
import tempfile
from typing import List, Optional

from chem_spectra.lib.composer.lcms_converter_app import LCMSConverterAppComposer
from chem_spectra.lib.external.chemotion_converter_lcms import (
    lcms_preview_image_from_jdx_files,
)

logger = logging.getLogger(__name__)


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


def build_lcms_composer(jcamp_paths, params=None) -> Optional[LCMSConverterAppComposer]:
    """Build a single LCMSConverterAppComposer from a group of LC/MS JCAMP files.

    The LC/MS pipeline expects ``NamedTemporaryFile`` handles, so each input
    path is copied into a fresh temporary file. Ownership of those handles is
    transferred to the returned composer (via ``composer.data``); the caller
    is responsible for closing them once consumed, which mirrors how
    ``BagItBaseConverter`` already manages NI/MS temp files (closed by the zip
    writer in ``to_zip_bag_it_response``).

    Returns ``None`` if the inputs are empty or if temp file creation fails so
    the caller can apply a safe fallback (skip the group).
    """
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
    """Build the LCMS composer for ``lcms_paths`` and extend the BagIt lists.

    Lists are mutated in place. On any failure the group is skipped so the
    rest of the BagIt dataset still produces a valid response.
    """
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

    if archive_stems is not None:
        ordered = sorted(lcms_paths)
        for path in ordered:
            base = os.path.basename(path)
            stem = os.path.splitext(base)[0].replace('.', '_')
            archive_stems.append(stem)

    list_files.extend(files)
    list_images.extend([preview] + [None] * (n - 1))
    list_csv.extend([csv_file] + [None] * (n - 1))
    list_composer.extend([composer] * n)
