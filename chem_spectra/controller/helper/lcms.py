"""LCMS-specific request helpers.

Keeps LCMS request parsing/normalisation out of the generic ``share.py``
helpers so the original (pre-LCMS) request pipeline stays untouched.
"""

import os.path as os_path


def _param_to_text(value):
    """Normalise an arbitrary form/file param into a UTF-8 string (or None)."""
    if value is None:
        return None

    if hasattr(value, 'read'):
        try:
            content = value.read()
        except Exception:
            return None
        try:
            value.seek(0)
        except Exception:
            pass
        if isinstance(content, bytes):
            return content.decode('utf-8', errors='ignore')
        return str(content)

    if isinstance(value, dict):
        tmp_file = value.get('tempfile') or value.get(':tempfile')
        if hasattr(tmp_file, 'read'):
            try:
                content = tmp_file.read()
            except Exception:
                return None
            try:
                tmp_file.seek(0)
            except Exception:
                pass
            if isinstance(content, bytes):
                return content.decode('utf-8', errors='ignore')
            return str(content)

    return value if isinstance(value, str) else str(value)


def normalize_lcms_filename(filename, src_filename=None):
    """Normalize an LCMS attachment filename for ELN grouping.

    When the frontend posts both a source archive (``src_filename``) and a
    derived ``filename`` (e.g. ``run_42.edit``), prefer the source basename so
    files from the same dataset stay grouped.

    Otherwise strip archive extensions and ``.edit`` suffixes from
    ``filename``. Intended for LCMS save paths only — callers should not
    invoke this helper for unrelated upload types.
    """
    if not filename and not src_filename:
        return filename

    def _strip_archive_suffix_simple(name):
        lower = name.lower()
        for suffix in (".tar.gz", ".tgz", ".tar.xz", ".tar", ".zip"):
            if lower.endswith(suffix):
                return name[: -len(suffix)]
        return os_path.splitext(name)[0]

    def _basename_no_ext(value):
        if not value:
            return None
        base = os_path.basename(str(value))
        return _strip_archive_suffix_simple(base)

    def _strip_edit_suffix(value):
        if not value:
            return value
        return value[:-5] if value.lower().endswith(".edit") else value

    base = _basename_no_ext(filename)
    src_base = _basename_no_ext(src_filename)

    if src_base:
        if not base:
            return src_base
        src_core = _strip_edit_suffix(src_base)
        base_core = _strip_edit_suffix(base)
        if base_core.lower().startswith(src_core.lower()):
            return src_base

    return base or filename


LCMS_PARAM_KEYS = (
    'lcms_uvvis_wavelength',
    'lcms_mz_page',
    'lcms_mz_page_data',
)


def extract_lcms_params(request):
    """Pull LCMS-specific form/file fields out of a Flask request.

    Returns a flat dict callers can merge with their general params. All keys
    default to ``None`` so non-LCMS requests stay untouched.
    """
    lcms_uvvis_wavelength = request.form.get('lcms_uvvis_wavelength', default=None)
    lcms_mz_page = request.form.get('lcms_mz_page', default=None)

    raw_mz_page_data = (
        request.files.get('lcms_mz_page_data')
        or request.form.get('lcms_mz_page_data', default=None)
    )
    lcms_mz_page_data = _param_to_text(raw_mz_page_data)

    return {
        'lcms_uvvis_wavelength': lcms_uvvis_wavelength,
        'lcms_mz_page': lcms_mz_page,
        'lcms_mz_page_data': lcms_mz_page_data,
    }
