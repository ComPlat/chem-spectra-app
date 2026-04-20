import io
import json
import os
import re
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


def _strip_archive_suffix(name: str) -> str:
    lower = name.lower()
    for suffix in (".tar.gz", ".tgz", ".tar.xz", ".tar", ".zip"):
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def _strip_lcms_tmp_prefix(name: str) -> str:
    if not name:
        return name
    lowered = name.lower()
    if lowered.startswith("tmp") and "_" in name:
        candidate = name.split("_", 1)[1]
        cand_lower = candidate.lower()
        if any(token in cand_lower for token in ("uvvis", "tic", "mz", "lcms", "chemstation")):
            return candidate
    return name


def _normalize_lcms_preview_base(name: str) -> Optional[str]:
    if not name:
        return None
    raw = os.path.basename(str(name))
    raw = _strip_lcms_tmp_prefix(raw)
    base = _strip_archive_suffix(raw)
    base = base.rstrip("_-")
    return base or None


def _lcms_preview_basename(paths: List[str], params: Optional[Dict]) -> Optional[str]:
    candidates: List[str] = []
    if isinstance(params, dict):
        for key in ("fname", "filename", "name"):
            value = params.get(key)
            if value:
                candidates.append(str(value))
    for entry in paths:
        candidates.append(entry)
    if candidates:
        normalized = [b for b in (_normalize_lcms_preview_base(c) for c in candidates) if b]
        if normalized:
            if len(normalized) > 1:
                prefix = os.path.commonprefix(normalized)
                if prefix:
                    return prefix.rstrip("_-")
            return normalized[0]
    return None


class _NamedFile:
    def __init__(self, fp, tmp_dir: tempfile.TemporaryDirectory):
        self._fp = fp
        self._tmp_dir = tmp_dir
        self.name = fp.name

    def __getattr__(self, name):
        return getattr(self._fp, name)

    def close(self):
        try:
            self._fp.close()
        finally:
            try:
                self._tmp_dir.cleanup()
            except Exception:
                pass


def _write_named_file(content: bytes, filename: str) -> _NamedFile:
    tmp_dir = tempfile.TemporaryDirectory(prefix="chemotion_lcms_")
    path = os.path.join(tmp_dir.name, filename)
    fp = open(path, "w+b")
    fp.write(content)
    fp.seek(0)
    return _NamedFile(fp, tmp_dir)


def _as_float(value) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


_FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")


def _float_from_label(value) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip().replace(",", ".")
    direct = _as_float(text)
    if direct is not None:
        return direct
    match = _FLOAT_RE.search(text)
    if not match:
        return None
    return _as_float(match.group(0))


def _format_number(value: Optional[float]) -> Optional[str]:
    if value is None:
        return None
    try:
        return f"{float(value):g}"
    except (TypeError, ValueError):
        return None


def _format_wavelength_label(value) -> Optional[str]:
    numeric = _float_from_label(value)
    if numeric is not None:
        return f"{_format_number(numeric)} nm"
    if value is None:
        return None
    text = str(value).strip()
    return f"{text} nm" if text else None


def _format_ms_page_label(value) -> Optional[str]:
    numeric = _float_from_label(value)
    if numeric is not None:
        return _format_number(numeric)
    if value is None:
        return None
    text = str(value).strip()
    return text if text else None


def _format_rt_label(value) -> Optional[str]:
    numeric = _float_from_label(value)
    if numeric is None:
        return None
    return f"{_format_number(numeric)} min"


def _extract_ms_peaks_from_param(value) -> List[Tuple[float, float]]:
    if not value:
        return []

    data = value
    if isinstance(data, str):
        raw = data.strip()
        if not raw:
            return []
        try:
            data = json.loads(raw)
        except Exception:
            return []

    if not isinstance(data, list):
        return []

    peaks: List[Tuple[float, float]] = []
    for item in data:
        if not isinstance(item, dict):
            continue
        mz = _as_float(item.get("x"))
        intensity = _as_float(item.get("y"))
        if mz is None or intensity is None:
            continue
        peaks.append((mz, intensity))
    return peaks


def _extract_xy_from_lines(lines: List[str], start_idx: int) -> Tuple[List[float], List[float]]:
    xs: List[float] = []
    ys: List[float] = []
    for line in lines[start_idx:]:
        stripped = line.strip()
        if not stripped:
            continue
        upper = stripped.upper()
        if upper.startswith("##") or upper.startswith("$$"):
            break
        numbers = _FLOAT_RE.findall(stripped.replace(";", " "))
        if len(numbers) < 2:
            continue
        for idx in range(0, len(numbers) - 1, 2):
            x_val = _as_float(numbers[idx])
            y_val = _as_float(numbers[idx + 1])
            if x_val is None or y_val is None:
                continue
            xs.append(x_val)
            ys.append(y_val)
    return xs, ys


def _extract_xy_from_jdx_content(content: str) -> Tuple[List[float], List[float]]:
    lines = content.splitlines()
    data_start = None
    for idx, line in enumerate(lines):
        upper = line.upper()
        if "##DATA TABLE" in upper or "##PEAK TABLE" in upper or "##XYDATA" in upper:
            data_start = idx + 1
            break
    if data_start is None:
        return [], []
    return _extract_xy_from_lines(lines, data_start)


def _extract_ms_page(
    content: str,
    target_page=None,
    *,
    first_page_only: bool = False,
) -> Tuple[List[float], List[float], Optional[float], Optional[str]]:
    """Parse MS JCAMP. If ``first_page_only`` and ``target_page`` is None, return after the first spectrum block."""
    lines = content.splitlines()
    threshold = None
    for line in lines:
        if line.upper().startswith("##$CSTHRESHOLD"):
            _, _, raw = line.partition("=")
            threshold = _as_float(raw.strip())
            break
    if threshold is not None and threshold > 1.0:
        threshold = threshold / 100.0
    pages: List[Tuple[str, List[float], List[float]]] = []
    current_label = None
    for idx, line in enumerate(lines):
        stripped = line.strip()
        upper = stripped.upper()
        if upper.startswith("##PAGE="):
            current_label = stripped.split("=", 1)[1].strip()
            continue
        if current_label and ("##DATA TABLE" in upper or "##PEAK TABLE" in upper or "##XYDATA" in upper):
            xs, ys = _extract_xy_from_lines(lines, idx + 1)
            if xs and ys:
                if first_page_only and target_page is None:
                    return xs, ys, threshold, current_label
                pages.append((current_label, xs, ys))
            continue

    if not pages:
        xs, ys = _extract_xy_from_jdx_content(content)
        return xs, ys, threshold, None

    target_label = None
    if target_page is not None:
        target_val = _float_from_label(target_page)
        if target_val is not None:
            numeric_pages = [(label, _float_from_label(label), xs, ys) for label, xs, ys in pages]
            numeric_only = [p for p in numeric_pages if p[1] is not None]
            if numeric_only:
                target_label = min(numeric_only, key=lambda p: abs(p[1] - target_val))[0]
        if target_label is None:
            target_str = str(target_page).strip()
            for label, _xs, _ys in pages:
                if str(label).strip() == target_str:
                    target_label = label
                    break

    for label, xs, ys in pages:
        if target_label is None or label == target_label:
            return xs, ys, threshold, label
    return pages[0][1], pages[0][2], threshold, pages[0][0]


def _mz_page_param_explicit(mz_page) -> bool:
    """True when the client asked for a specific MS page / RT (not default first page)."""
    if mz_page is None:
        return False
    return bool(str(mz_page).strip())


def _extract_ms_page_first_page_from_path(path: str) -> Tuple[List[float], List[float], Optional[float], Optional[str]]:
    """Read MZ JCAMP line-by-line; stop after the first spectrum block (preview / fresh upload).

    Avoids loading multi-hundred-MB files and parsing every ##PAGE block.
    """
    threshold: Optional[float] = None
    lines: List[str] = []
    current_label: Optional[str] = None
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as handle:
            for raw in handle:
                line = raw.rstrip("\n\r")
                lines.append(line)
                idx = len(lines) - 1
                lu = line.upper()
                if lu.startswith("##$CSTHRESHOLD"):
                    _, _, r = line.partition("=")
                    threshold = _as_float(r.strip())
                    continue
                stripped = line.strip()
                upper = stripped.upper()
                if upper.startswith("##PAGE="):
                    current_label = stripped.split("=", 1)[1].strip()
                    continue
                if current_label and (
                    "##DATA TABLE" in upper
                    or "##PEAK TABLE" in upper
                    or "##XYDATA" in upper
                ):
                    xs, ys = _extract_xy_from_lines(lines, idx + 1)
                    if xs and ys:
                        thr = threshold
                        if thr is not None and thr > 1.0:
                            thr = thr / 100.0
                        return xs, ys, thr, current_label
    except OSError:
        pass

    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as handle:
            content = handle.read()
    except OSError:
        return [], [], None, None
    return _extract_ms_page(content, None, first_page_only=True)


def _extract_uvvis_from_peak_content(content: str, target_wavelength=None):
    if "CHEMSPECTRA UVVIS PEAK TABLE" not in content:
        return None
    try:
        import numpy as np  # type: ignore
    except Exception:
        return None

    data_by_wavelength = {}
    peaks_by_wavelength = {}
    integrations_by_wavelength = {}

    sections = content.split("$$ === CHEMSPECTRA UVVIS PEAK TABLE ===")
    for section in sections[1:]:
        wavelength = None
        for line in section.split("\n"):
            if line.startswith("##PAGE="):
                try:
                    wavelength = float(line.replace("##PAGE=", "").strip())
                except (ValueError, AttributeError):
                    wavelength = line.replace("##PAGE=", "").strip()
                break
        if wavelength is None:
            continue

        xs = []
        ys = []
        in_data_section = False
        for line in section.split("\n"):
            if "##DATA TABLE=" in line.upper():
                in_data_section = True
                continue
            if "##END=" in line or "$$" in line:
                break
            if in_data_section and line.strip() and not line.strip().startswith("##"):
                cleaned = line.replace(";", "").strip()
                if not cleaned:
                    continue
                parts = cleaned.split(",")
                if len(parts) >= 2:
                    x_val = _as_float(parts[0].strip())
                    y_val = _as_float(parts[1].strip())
                    if x_val is not None and y_val is not None:
                        xs.append(x_val)
                        ys.append(y_val)

        if xs and ys:
            data_by_wavelength[wavelength] = (xs, ys)

    edit_sections = content.split("$$ === CHEMSPECTRA PEAK TABLE EDIT ===")
    for section in edit_sections[1:]:
        wavelength = None
        for line in section.split("\n"):
            if line.startswith("##PAGE="):
                try:
                    wavelength = float(line.replace("##PAGE=", "").strip())
                except (ValueError, AttributeError):
                    wavelength = line.replace("##PAGE=", "").strip()
                break
        if wavelength is None:
            continue

        peaks = []
        in_peak_section = False
        for line in section.split("\n"):
            if "##PEAKTABLE=" in line.upper():
                in_peak_section = True
                continue
            if "##END=" in line or "$$" in line:
                break
            if in_peak_section and line.strip() and not line.strip().startswith("##"):
                parts = line.strip().split(",")
                if len(parts) >= 2:
                    x_val = _as_float(parts[0].strip())
                    y_val = _as_float(parts[1].strip())
                    if x_val is not None and y_val is not None:
                        peaks.append({"x": x_val, "y": y_val})
        if peaks:
            peaks_by_wavelength[wavelength] = peaks

    intg_sections = content.split("$$ === CHEMSPECTRA INTEGRALS AND MULTIPLETS ===")
    for section in intg_sections[1:]:
        wavelength = None
        for line in section.split("\n"):
            if line.startswith("##PAGE="):
                try:
                    wavelength = float(line.replace("##PAGE=", "").strip())
                except (ValueError, AttributeError):
                    wavelength = line.replace("##PAGE=", "").strip()
                break
        if wavelength is None:
            continue

        integrations = []
        in_intg_section = False
        for line in section.split("\n"):
            if "##$OBSERVEDINTEGRALS=" in line.upper():
                in_intg_section = True
                continue
            if "##END=" in line or "$$" in line:
                break
            if in_intg_section and line.strip() and not line.strip().startswith("##"):
                cleaned = line.strip().replace("(", "").replace(")", "")
                parts = cleaned.split(",")
                if len(parts) >= 3:
                    x_left = _as_float(parts[0].strip())
                    x_right = _as_float(parts[1].strip())
                    area = _as_float(parts[2].strip())
                    if x_left is not None and x_right is not None and area is not None:
                        integrations.append({"xL": x_left, "xU": x_right, "area": area})
        if integrations:
            integrations_by_wavelength[wavelength] = integrations

    if not data_by_wavelength:
        return None

    def _to_float_or_none(value):
        try:
            return float(value)
        except Exception:
            return None

    numeric_pairs = [(wl, _to_float_or_none(wl)) for wl in data_by_wavelength.keys()]
    numeric_only = [p for p in numeric_pairs if p[1] is not None]

    wl_key = None
    if target_wavelength is not None:
        if target_wavelength in data_by_wavelength:
            wl_key = target_wavelength
        else:
            target_val = _to_float_or_none(target_wavelength)
            if target_val is not None and numeric_only:
                wl_key = min(numeric_only, key=lambda p: abs(p[1] - target_val))[0]
    if wl_key is None:
        wl_key = min(numeric_only, key=lambda p: p[1])[0] if numeric_only else list(data_by_wavelength.keys())[0]

    xs, ys = data_by_wavelength[wl_key]
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)

    edit_peaks = peaks_by_wavelength.get(wl_key, [])
    integrations = integrations_by_wavelength.get(wl_key, [])
    return xs, ys, edit_peaks, integrations, wl_key


def _pick_jdx_path(paths: List[str], primary_token: str, fallback_token: Optional[str] = None) -> Optional[str]:
    lowered = [(p, os.path.basename(p).lower()) for p in paths]
    for p, name in lowered:
        if primary_token in name:
            return p
    if fallback_token:
        for p, name in lowered:
            if fallback_token in name:
                return p
    return None


def _header_from_content(content: str, limit: int = 200) -> Dict[str, str]:
    header: Dict[str, str] = {}
    for line in content.splitlines()[:limit]:
        if not line.startswith("##"):
            continue
        if "=" not in line:
            continue
        key, value = line[2:].split("=", 1)
        header[key.strip()] = value.strip()
    return header


def _classify_lcms_content(content: str) -> Optional[str]:
    header = _header_from_content(content)
    kind, _ = _classify_lcms_header(header)
    return kind


def _normalize_params(params: Optional[Dict]) -> Dict:
    if not params:
        return {}
    if isinstance(params, dict):
        integration = params.get("integration")
        if isinstance(integration, dict):
            return params
    try:
        from chem_spectra.lib.converter.share import parse_params  # type: ignore
    except Exception:
        return params if isinstance(params, dict) else {}
    try:
        return parse_params(params)
    except Exception:
        return params if isinstance(params, dict) else {}


def _extract_uvvis_pages(content: str) -> List[Tuple[Optional[float], List[float], List[float]]]:
    lines = content.splitlines()
    pages: List[Tuple[Optional[float], List[float], List[float]]] = []
    current_page: Optional[float] = None
    for idx, line in enumerate(lines):
        stripped = line.strip()
        upper = stripped.upper()
        if upper.startswith("##PAGE="):
            current_page = _float_from_label(stripped.split("=", 1)[1].strip())
            continue
        if "##DATA TABLE" in upper or "##PEAK TABLE" in upper or "##XYDATA" in upper:
            xs, ys = _extract_xy_from_lines(lines, idx + 1)
            if xs and ys:
                pages.append((current_page, xs, ys))
    if pages:
        return pages

    xs, ys = _extract_xy_from_jdx_content(content)
    if not xs or not ys:
        return []
    return [(None, xs, ys)]


def _lc_df_from_uvvis_content(content: str) -> Optional[pd.DataFrame]:
    pages = _extract_uvvis_pages(content)
    if not pages:
        return None
    rows = []
    for wavelength, xs, ys in pages:
        page_wavelength = wavelength if wavelength is not None else 0.0
        for x_val, y_val in zip(xs, ys):
            rows.append(
                {
                    "RetentionTime": x_val,
                    "DetectorSignal": y_val,
                    "wavelength": page_wavelength,
                }
            )
    if not rows:
        return None
    return pd.DataFrame(rows, columns=["RetentionTime", "DetectorSignal", "wavelength"])


def _header_value(header: Dict, key: str) -> Optional[str]:
    for k, v in header.items():
        if str(k).strip().lower() == key.lower():
            return str(v)
    return None


def _classify_lcms_header(header: Dict) -> Tuple[Optional[str], Optional[str]]:
    category = _header_value(header, "$CSCATEGORY") or _header_value(header, "CSCATEGORY")
    ntuples_id = _header_value(header, "NTUPLES_ID")
    data_type = _header_value(header, "DATA TYPE")
    data_class = _header_value(header, "DATA CLASS")
    x_units = _header_value(header, "XUNITS")
    y_units = _header_value(header, "YUNITS")
    scan_mode = _header_value(header, "SCAN_MODE") or _header_value(header, "$SCAN_MODE")
    ion_mode = _header_value(header, "ION_MODE") or _header_value(header, "$ION_MODE")
    mode = _header_value(header, "MODE") or _header_value(header, "$MODE")
    polarity_header = _header_value(header, "POLARITY") or _header_value(header, "$POLARITY")
    combined = " ".join([v for v in [category, ntuples_id, data_type] if v])
    combined_upper = combined.upper()

    kind: Optional[str] = None
    if "UVVIS" in combined_upper or "UV/VIS" in combined_upper or "UV-VIS" in combined_upper or "UV VIS" in combined_upper:
        kind = "uvvis"
    elif "TIC" in combined_upper:
        kind = "tic"
    elif "MZ" in combined_upper or "M/Z" in combined_upper or "MASS SPECTRUM" in combined_upper:
        kind = "mz"
    if not kind and data_class and "XYPOINTS" in data_class.upper():
        data_type_upper = data_type.upper() if data_type else ""
        x_units_upper = x_units.upper() if x_units else ""
        y_units_upper = y_units.upper() if y_units else ""
        if "LC/MS" in data_type_upper or "MASS TIC" in data_type_upper or "TIC" in data_type_upper:
            kind = "tic"
        elif "MINUTE" in x_units_upper and "INTENS" in y_units_upper:
            kind = "tic"
    if not kind:
        return None, None

    polarity_combined_upper = " ".join(
        [v for v in [combined, scan_mode, ion_mode, mode, polarity_header] if v]
    ).upper()
    polarity: Optional[str] = None
    if "NEG" in polarity_combined_upper or "MINUS" in polarity_combined_upper:
        polarity = "minus"
    elif "POS" in polarity_combined_upper or "PLUS" in polarity_combined_upper:
        polarity = "plus"
    return kind, polarity


def lcms_preview_image_from_jdx_files(
    jdx_files: List,
    params: Optional[Dict] = None,
) -> Optional[tempfile.NamedTemporaryFile]:
    """Build LC/MS preview PNG."""
    if not jdx_files:
        return None
    paths: List[str] = []
    for entry in jdx_files:
        if isinstance(entry, str):
            paths.append(entry)
        else:
            path = getattr(entry, "name", None)
            if path:
                paths.append(path)
    if not paths:
        return None

    normalized_params = _normalize_params(params)
    uvvis_wavelength = normalized_params.get("lcms_uvvis_wavelength")
    mz_page = normalized_params.get("lcms_mz_page")
    mz_page_data = normalized_params.get("lcms_mz_page_data")
    ms_threshold_from_param = _as_float(normalized_params.get("thres"))
    selected_ms_peaks = _extract_ms_peaks_from_param(mz_page_data)

    peak_path = _pick_jdx_path(paths, "uvvis.peak", "peak.jdx")
    mz_path = _pick_jdx_path(paths, "_mz", "mz")

    if not peak_path:
        for path in paths:
            try:
                with open(path, "r", encoding="utf-8", errors="ignore") as handle:
                    content = handle.read()
                if "CHEMSPECTRA UVVIS PEAK TABLE" in content:
                    peak_path = path
                    break
            except Exception:
                continue

    if not mz_path:
        for path in paths:
            if path == peak_path:
                continue
            try:
                with open(path, "r", encoding="utf-8", errors="ignore") as handle:
                    content = handle.read()
                kind = _classify_lcms_content(content)
                if kind == "mz" and not mz_path:
                    mz_path = path
            except Exception:
                continue

    uvvis_data = None
    uvvis_label = None
    if peak_path:
        try:
            with open(peak_path, "r", encoding="utf-8", errors="ignore") as handle:
                uvvis_content = handle.read()
            uvvis_data = _extract_uvvis_from_peak_content(uvvis_content, uvvis_wavelength)
            if uvvis_data:
                uvvis_label = _format_wavelength_label(uvvis_data[4])
        except Exception:
            uvvis_data = None

    ms_data = None
    ms_threshold = None
    ms_page_label = _format_ms_page_label(mz_page)
    ms_rt = _float_from_label(mz_page)
    ms_rt_label = _format_rt_label(ms_rt)
    if mz_path:
        try:
            if selected_ms_peaks:
                # Pics MS déjà fournis par l'ELN : inutile de lire tout le JCAMP MZ.
                ms_data = None
            elif _mz_page_param_explicit(mz_page):
                with open(mz_path, "r", encoding="utf-8", errors="ignore") as handle:
                    mz_content = handle.read()
                ms_data = _extract_ms_page(mz_content, mz_page, first_page_only=False)
            else:
                # Nouvel upload : uniquement la première page, lecture flux (arrêt après 1er spectre).
                ms_data = _extract_ms_page_first_page_from_path(mz_path)

            if ms_data and ms_data[0]:
                ms_threshold = ms_data[2]
                ms_page_label = _format_ms_page_label(ms_data[3]) or _format_ms_page_label(mz_page)
                if ms_rt is None:
                    ms_rt = _float_from_label(ms_data[3])
                    ms_rt_label = _format_rt_label(ms_rt)
        except Exception:
            ms_data = None

    active_ms_threshold = ms_threshold_from_param
    if active_ms_threshold is None or active_ms_threshold <= 0.0:
        active_ms_threshold = ms_threshold
    if active_ms_threshold is not None and active_ms_threshold > 1.0:
        active_ms_threshold = active_ms_threshold / 100.0
    if active_ms_threshold is None or active_ms_threshold <= 0.0:
        # Default to 5% when threshold is missing/uninitialized.
        active_ms_threshold = 0.05
    active_ms_threshold = max(0.0, min(1.0, active_ms_threshold))
    has_uvvis = uvvis_data is not None
    has_ms = (ms_data is not None and ms_data[0]) or bool(selected_ms_peaks)
    if not (has_uvvis or has_ms):
        return None

    try:
        import matplotlib.pyplot as plt  # type: ignore
        import matplotlib.path as mpath  # type: ignore
        from matplotlib.ticker import MultipleLocator  # type: ignore
        import numpy as np  # type: ignore
        from chem_spectra.lib.shared.calc import get_curve_endpoint, cal_slope
    except Exception:
        return None

    plt.rcParams["figure.figsize"] = [16, 9]
    plt.rcParams["figure.dpi"] = 200
    plt.rcParams["font.size"] = 14

    fig, axes = plt.subplots(2, 1)
    uvvis_ax, ms_ax = axes

    if has_uvvis:
        xs, ys, edit_peaks, integrations, _wl = uvvis_data
        xs_min = xs / 60.0
        uvvis_ax.plot(xs_min, ys)

        if edit_peaks:
            path_data = [
                (mpath.Path.MOVETO, (0, 5)),
                (mpath.Path.LINETO, (0, 20)),
            ]
            codes, verts = zip(*path_data)
            marker = mpath.Path(verts, codes)
            x_peaks = [p["x"] / 60.0 for p in edit_peaks]
            y_peaks = [p["y"] for p in edit_peaks]
            uvvis_ax.plot(
                x_peaks,
                y_peaks,
                "r",
                ls="",
                marker=marker,
                markersize=50,
            )

        if integrations:
            y_max = float(np.max(ys))
            y_min = float(np.min(ys))
            h = max(y_max - y_min, 1.0)
            for itg in integrations:
                x_left = itg["xL"] / 60.0
                x_right = itg["xU"] / 60.0
                i_left, i_right = get_curve_endpoint(xs_min, ys, x_left, x_right)
                cxs = xs_min[i_left:i_right]
                cys = ys[i_left:i_right]
                if len(cxs) > 0 and len(cys) > 0:
                    slope = cal_slope(cxs[0], cys[0], cxs[len(cxs) - 1], cys[len(cys) - 1])
                    last_y = cys[0]
                    last_x = cxs[0]
                    aucys = [last_y]
                    for i in range(1, len(cys)):
                        curr_x = cxs[i]
                        curr_y = slope * (curr_x - last_x) + last_y
                        aucys.append(curr_y)
                        last_x = curr_x
                        last_y = curr_y
                    uvvis_ax.fill_between(cxs, y1=cys, y2=aucys, alpha=0.2, color="#FF0000")

        uvvis_ax.set_xlabel("X (Retention Time, min)", fontsize=18)
        uvvis_ax.set_ylabel("Y (Detector Signal)", fontsize=18)
        uvvis_ax.grid(False)
        if xs_min.size:
            x_min, x_max = float(np.min(xs_min)), float(np.max(xs_min))
        else:
            x_min, x_max = 0.0, 1.0
        if ys.size:
            y_min, y_max = float(np.min(ys)), float(np.max(ys))
        else:
            y_min, y_max = 0.0, 1.0
        h = max(y_max - y_min, 1.0)
        uvvis_ax.set_xlim(x_min, x_max)
        uvvis_ax.xaxis.set_major_locator(MultipleLocator(1.0))
        uvvis_ax.set_ylim(y_min - h * 0.05, y_max + h * 0.15)
        if uvvis_label:
            uvvis_ax.text(
                0.98,
                0.98,
                uvvis_label,
                ha="right",
                va="top",
                transform=uvvis_ax.transAxes,
            )
        if ms_rt is not None:
            uvvis_ax.axvline(ms_rt, color="#777777", linestyle="--", linewidth=1, alpha=0.7)
            if ms_rt_label:
                uvvis_ax.text(
                    0.02,
                    0.98,
                    f"MS RT {ms_rt_label}",
                    ha="left",
                    va="top",
                    transform=uvvis_ax.transAxes,
                    fontsize=12,
                    color="#555555",
                )
    else:
        uvvis_ax.text(0.5, 0.5, "UVVIS unavailable", ha="center", va="center", transform=uvvis_ax.transAxes)
        uvvis_ax.set_axis_off()

    if has_ms:
        ms_panel_xs: List[float] = []
        if ms_data is not None and ms_data[0]:
            ms_xs, ms_ys, _ms_threshold, _page = ms_data
            if active_ms_threshold is not None and ms_ys:
                cut = max(ms_ys) * active_ms_threshold
                ms_xs_blue = []
                ms_ys_blue = []
                ms_xs_grey = []
                ms_ys_grey = []
                for ms_x, ms_y in zip(ms_xs, ms_ys):
                    if ms_y >= cut:
                        ms_xs_blue.append(ms_x)
                        ms_ys_blue.append(ms_y)
                    else:
                        ms_xs_grey.append(ms_x)
                        ms_ys_grey.append(ms_y)
                if ms_xs_grey:
                    ms_ax.bar(ms_xs_grey, ms_ys_grey, width=0, edgecolor="#dddddd")
                if ms_xs_blue:
                    ms_ax.bar(ms_xs_blue, ms_ys_blue, width=0, edgecolor="#1f77b4")
            else:
                ms_ax.bar(ms_xs, ms_ys, width=0, edgecolor="#dddddd")
            ms_panel_xs.extend(ms_xs)

        if selected_ms_peaks:
            peak_xs = [mz for mz, _ in selected_ms_peaks]
            peak_ys = [intensity for _, intensity in selected_ms_peaks]
            if active_ms_threshold is not None and peak_ys:
                cut = max(peak_ys) * active_ms_threshold
                peak_xs_blue = []
                peak_ys_blue = []
                peak_xs_grey = []
                peak_ys_grey = []
                for peak_x, peak_y in zip(peak_xs, peak_ys):
                    if peak_y >= cut:
                        peak_xs_blue.append(peak_x)
                        peak_ys_blue.append(peak_y)
                    else:
                        peak_xs_grey.append(peak_x)
                        peak_ys_grey.append(peak_y)
                if peak_xs_grey:
                    ms_ax.bar(peak_xs_grey, peak_ys_grey, width=0, edgecolor="#dddddd")
                if peak_xs_blue:
                    ms_ax.bar(peak_xs_blue, peak_ys_blue, width=0, edgecolor="#1f77b4")
            else:
                ms_ax.bar(peak_xs, peak_ys, width=0, edgecolor="#1f77b4")
            ms_panel_xs.extend(peak_xs)
        ms_ax.set_xlabel("X (m/z)", fontsize=18)
        ms_ax.set_ylabel("Y (Relative Abundance)", fontsize=18)
        ms_ax.grid(False)
        if ms_panel_xs:
            ms_ax.set_xlim(min(ms_panel_xs), max(ms_panel_xs))
            ms_ax.margins(x=0)
        if ms_rt_label:
            ms_ax.text(
                0.98,
                0.98,
                f"RT {ms_rt_label}",
                ha="right",
                va="top",
                transform=ms_ax.transAxes,
            )
        elif ms_page_label:
            ms_ax.text(
                0.98,
                0.98,
                f"page {ms_page_label}",
                ha="right",
                va="top",
                transform=ms_ax.transAxes,
            )
    else:
        ms_ax.text(0.5, 0.5, "MS unavailable", ha="center", va="center", transform=ms_ax.transAxes)
        ms_ax.set_axis_off()

    fig.tight_layout()
    preview_base = _lcms_preview_basename(paths, params)
    preview_name = f"{preview_base}.png" if preview_base else "lcms_preview.png"
    buffer = io.BytesIO()
    plt.savefig(buffer, format="png")
    buffer.seek(0)
    tf_img = _write_named_file(buffer.read(), preview_name)
    buffer.close()
    plt.clf()
    plt.cla()
    plt.close(fig)
    return tf_img


def lcms_df_from_peak_jdx(jdx_path: str) -> Optional[pd.DataFrame]:
    import pandas as pd
    
    if not jdx_path or not os.path.exists(jdx_path):
        return None
    
    try:
        with open(jdx_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        
        if 'CHEMSPECTRA UVVIS PEAK TABLE' not in content:
            return None
        
        rows = []
        sections = content.split('$$ === CHEMSPECTRA UVVIS PEAK TABLE ===')
        for section in sections[1:]:
            wavelength = None
            for line in section.split('\n'):
                if line.startswith('##PAGE='):
                    try:
                        wavelength = float(line.replace('##PAGE=', '').strip())
                    except (ValueError, AttributeError):
                        wavelength = line.replace('##PAGE=', '').strip()
                    break
            
            if wavelength is None:
                continue
            
            xs = []
            ys = []
            in_data_section = False
            for line in section.split('\n'):
                if '##DATA TABLE=' in line.upper():
                    in_data_section = True
                    continue
                if '##END=' in line or '$$' in line:
                    break
                if in_data_section and line.strip() and not line.strip().startswith('##'):
                    cleaned = line.replace(';', '').strip()
                    if not cleaned:
                        continue
                    parts = cleaned.split(',')
                    if len(parts) >= 2:
                        try:
                            x = _as_float(parts[0].strip())
                            y = _as_float(parts[1].strip())
                            if x is not None and y is not None:
                                xs.append(x)
                                ys.append(y)
                        except (ValueError, IndexError):
                            continue
            
            for x, y in zip(xs, ys):
                rows.append({
                    'RetentionTime': x,
                    'DetectorSignal': y,
                    'wavelength': wavelength,
                })
        
        if rows:
            df = pd.DataFrame(rows, columns=['RetentionTime', 'DetectorSignal', 'wavelength'])
            return df
        return None
    except Exception as e:
        return None


def lcms_df_from_uvvis_jdx(jdx_path: str) -> Optional[pd.DataFrame]:
    if not jdx_path or not os.path.exists(jdx_path):
        return None
    try:
        with open(jdx_path, "r", encoding="utf-8", errors="ignore") as handle:
            content = handle.read()
        return _lc_df_from_uvvis_content(content)
    except Exception:
        return None

def lcms_uvvis_peak_jcamp_from_df(
    lc_df: pd.DataFrame,
    title: str,
    params: Optional[Dict] = None,
) -> Optional[tempfile.NamedTemporaryFile]:
    if lc_df is None or lc_df.empty or "wavelength" not in lc_df.columns:
        return None

    normalized_params = _normalize_params(params)
    peaks_input = normalized_params.get("peaks_str") if normalized_params else None
    integration_input = normalized_params.get("integration") if normalized_params else None
    selected_page_idx = normalized_params.get("jcamp_idx", 0) if normalized_params else 0
    lcms_mz_page_hdr = ""
    if normalized_params:
        _raw_mz_page = normalized_params.get("lcms_mz_page")
        if _raw_mz_page is not None:
            lcms_mz_page_hdr = " ".join(str(_raw_mz_page).splitlines()).strip()

    def _to_float_or_none(val):
        try:
            return float(val)
        except Exception:
            return None

    def _normalize_peaks_list(value) -> List[Dict[str, float]]:
        if not value:
            return []
        if isinstance(value, dict):
            if "x" in value and "y" in value:
                xs = value.get("x") or []
                ys = value.get("y") or []
                return [{"x": x, "y": y} for x, y in zip(xs, ys)]
            return []
        if isinstance(value, list):
            if all(isinstance(v, dict) and "x" in v and "y" in v for v in value):
                return [{"x": v["x"], "y": v["y"]} for v in value]
            peaks = []
            for item in value:
                if isinstance(item, (list, tuple)) and len(item) >= 2:
                    peaks.append({"x": item[0], "y": item[1]})
            return peaks
        if isinstance(value, str):
            items = []
            for p in value.split("#"):
                parts = p.split(",")
                if len(parts) < 2:
                    continue
                x_val = _to_float_or_none(parts[0])
                y_val = _to_float_or_none(parts[1])
                if x_val is None or y_val is None:
                    continue
                items.append({"x": x_val, "y": y_val})
            return items
        return []

    def _parse_peaks_input(value):
        if value is None:
            return None
        if isinstance(value, (dict, list)):
            return value
        if isinstance(value, str):
            raw = value.strip()
            if raw.startswith("{") or raw.startswith("["):
                try:
                    return json.loads(raw)
                except Exception:
                    return value
        return value

    def _key_matches_page(key, page_idx, wavelength_value) -> bool:
        if key is None:
            return False
        if str(key) == str(page_idx):
            return True
        if str(key) == str(wavelength_value):
            return True
        key_num = _to_float_or_none(key)
        wl_num = _to_float_or_none(wavelength_value)
        return key_num is not None and wl_num is not None and key_num == wl_num

    def _peaks_for_page(peaks_data, page_idx, wavelength_value) -> List[Dict[str, float]]:
        if not peaks_data:
            return []
        peaks_data = _parse_peaks_input(peaks_data)
        if isinstance(peaks_data, dict):
            if "x" in peaks_data and "y" in peaks_data:
                return _normalize_peaks_list(peaks_data) if page_idx == selected_page_idx else []
            for key, value in peaks_data.items():
                if _key_matches_page(key, page_idx, wavelength_value):
                    return _normalize_peaks_list(value)
            return []
        if isinstance(peaks_data, list):
            if all(isinstance(v, dict) and "x" in v and "y" in v for v in peaks_data):
                return peaks_data if page_idx == selected_page_idx else []
            for item in peaks_data:
                if not isinstance(item, dict):
                    continue
                curve_idx = item.get("curveIdx") or item.get("curve_idx")
                curve_wl = item.get("wavelength")
                if curve_idx is not None and int(curve_idx) == int(page_idx):
                    return _normalize_peaks_list(item.get("peaks") or item.get("points") or item)
                if curve_wl is not None and _key_matches_page(curve_wl, page_idx, wavelength_value):
                    return _normalize_peaks_list(item.get("peaks") or item.get("points") or item)
            return []
        if isinstance(peaks_data, str):
            return _normalize_peaks_list(peaks_data) if page_idx == selected_page_idx else []
        return []

    def _integration_for_page(integration_data, page_idx, wavelength_value):
        if not integration_data:
            return None, []
        if isinstance(integration_data, list):
            for item in integration_data:
                if not isinstance(item, dict):
                    continue
                curve_idx = item.get("curveIdx") or item.get("curve_idx")
                curve_wl = item.get("wavelength")
                if curve_idx is not None and int(curve_idx) == int(page_idx):
                    return item, item.get("stack", []) or []
                if curve_wl is not None and _key_matches_page(curve_wl, page_idx, wavelength_value):
                    return item, item.get("stack", []) or []
            return None, []
        if isinstance(integration_data, dict):
            # Check if it's a dict with wavelength keys like {"230":[[xL, xU, area, absoluteArea]]}
            for key, value in integration_data.items():
                if _key_matches_page(key, page_idx, wavelength_value):
                    # Convert list format [[xL, xU, area, absoluteArea]] to dict format
                    stack = []
                    if isinstance(value, list):
                        for item in value:
                            if isinstance(item, (list, tuple)) and len(item) >= 3:
                                itg_dict = {
                                    'xL': item[0],
                                    'xU': item[1],
                                    'area': item[2],
                                    'absoluteArea': item[3] if len(item) >= 4 else 0
                                }
                                stack.append(itg_dict)
                            elif isinstance(item, dict):
                                stack.append(item)
                    return integration_data, stack
            
            curves = integration_data.get("curves") or integration_data.get("byCurve") or integration_data.get("by_curve")
            if isinstance(curves, list):
                for curve in curves:
                    if not isinstance(curve, dict):
                        continue
                    curve_idx = curve.get("curveIdx") or curve.get("curve_idx")
                    curve_wl = curve.get("wavelength")
                    if curve_idx is not None and int(curve_idx) == int(page_idx):
                        return curve, curve.get("stack", []) or []
                    if curve_wl is not None and _key_matches_page(curve_wl, page_idx, wavelength_value):
                        return curve, curve.get("stack", []) or []
                return None, []
            if page_idx == selected_page_idx:
                return integration_data, integration_data.get("stack", []) or []
        return None, []

    def _integration_lines(integration_meta, stack):
        if not stack:
            return []
        ref_area = 1.0
        ref_shift = 0.0
        if isinstance(integration_meta, dict):
            ref_area_raw = integration_meta.get("refArea") or 1
            ref_factor_raw = integration_meta.get("refFactor") or 1
            ref_shift = integration_meta.get("shift") or 0
            try:
                ref_area = float(ref_factor_raw) / float(ref_area_raw)
            except Exception:
                ref_area = 1.0
        lines = []
        for itg in stack:
            if not isinstance(itg, dict):
                continue
            x_left = itg.get("xL")
            x_right = itg.get("xU")
            area = itg.get("area")
            absolute_area = itg.get("absoluteArea", 0)
            if x_left is None or x_right is None or area is None:
                continue
            try:
                # Convert from minutes to seconds
                x_left = float(x_left) * 60.0 - float(ref_shift)
                x_right = float(x_right) * 60.0 - float(ref_shift)
                area = float(area) * ref_area
                absolute_area = float(absolute_area)
            except Exception:
                continue
            lines.append(f"({x_left}, {x_right}, {area}, {absolute_area})\n")
        return lines

    wavelength_values = lc_df["wavelength"].dropna().unique().tolist()
    numeric_pairs = [(w, _to_float_or_none(w)) for w in wavelength_values]
    numeric_only = [p for p in numeric_pairs if p[1] is not None]
    if numeric_only:
        ordered = [p[0] for p in sorted(numeric_only, key=lambda p: p[1])]
    else:
        ordered = sorted(wavelength_values)

    content: List[str] = []
    for page_idx, wl in enumerate(ordered):
        group = lc_df[lc_df["wavelength"] == wl]
        xs = group["RetentionTime"].tolist()
        ys = group["DetectorSignal"].tolist()
        if xs:
            x_min, x_max = min(xs), max(xs)
        else:
            x_min, x_max = 0.0, 0.0
        if ys:
            y_min, y_max = min(ys), max(ys)
        else:
            y_min, y_max = 0.0, 0.0

        edit_peaks = _peaks_for_page(peaks_input, page_idx, wl)
        # Convert peak x values from minutes to seconds
        for peak in edit_peaks:
            if 'x' in peak:
                peak['x'] = peak['x'] * 60.0
        
        integration_meta, integration_stack = _integration_for_page(integration_input, page_idx, wl)
        integration_lines = _integration_lines(integration_meta or integration_input, integration_stack)

        page_block = [
            "\n",
            "$$ === CHEMSPECTRA UVVIS PEAK TABLE ===\n",
            f"##TITLE={title}\n",
            "##JCAMP-DX=5.00\n",
            "##DATA TYPE=LC/MS\n",
            "##DATA CLASS=PEAK TABLE\n",
            "##SYMBOL=X, Y\n",
            "##ORIGIN=\n",
            "##OWNER=\n",
            "##XUNITS=RETENTION TIME\n",
            "##YUNITS=DETECTOR SIGNAL\n",
            "##$CSCATEGORY=UVVIS PEAK TABLE\n",
        ]
        if page_idx == 0 and lcms_mz_page_hdr:
            page_block.append(f"##$CSLCMSMZPAGE={lcms_mz_page_hdr}\n")
        page_block.extend([
            f"##PAGE={wl}\n",
            f"##NPOINTS={len(xs)}\n",
            "##DATA TABLE= (XY..XY), PEAKS\n",
        ])
        for x, y in zip(xs, ys):
            page_block.append(f"{x}, {y};\n")
        page_block.append("##END=\n")
        content.extend(page_block)

        content.extend([
            "\n",
            "$$ === CHEMSPECTRA INTEGRALS AND MULTIPLETS ===\n",
            f"##PAGE={wl}\n",
            "##$OBSERVEDINTEGRALS= (X Y Z)\n",
        ])
        content.extend(integration_lines)
        content.extend([
            "##$OBSERVEDMULTIPLETS=\n",
            "##$OBSERVEDMULTIPLETSPEAKS=\n",
            "##END=\n",
        ])

        content.extend([
            "\n",
            "$$ === CHEMSPECTRA PEAK TABLE EDIT ===\n",
            f"##TITLE={title}\n",
            "##JCAMP-DX=5.00\n",
            "##DATA TYPE=LC/MSPEAKTABLE\n",
            "##DATA CLASS=PEAKTABLE\n",
            "##$CSCATEGORY=EDIT_PEAK\n",
            "##$CSTHRESHOLD=0.05\n",
            f"##MAXX={x_max}\n",
            f"##MAXY={y_max}\n",
            f"##MINX={x_min}\n",
            f"##MINY={y_min}\n",
            f"##PAGE={wl}\n",
            f"##NPOINTS={len(edit_peaks)}\n",
            "##PEAKTABLE= (XY..XY)\n",
        ])
        for peak in edit_peaks:
            content.append(f"{peak['x']}, {peak['y']}\n")
        content.extend(["##END=\n"])

    tf = tempfile.NamedTemporaryFile(suffix="_lcms_uvvis.peak.jdx")
    tf.write(bytes("".join(content), "utf-8"))
    tf.seek(0)
    return tf

