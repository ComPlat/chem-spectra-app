import json
import os
import tarfile
import tempfile
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

def _find_first_cdf(root_dir: str) -> Optional[str]:
    for dirpath, _, filenames in os.walk(root_dir):
        for fn in filenames:
            if fn.lower().endswith(".cdf"):
                return os.path.join(dirpath, fn)
    return None


def _is_tar_path(path: str) -> bool:
    lower = path.lower()
    return lower.endswith(".tar") or lower.endswith(".tar.gz") or lower.endswith(".tgz") or lower.endswith(".tar.xz")


def _strip_archive_suffix(name: str) -> str:
    lower = name.lower()
    for suffix in (".tar.gz", ".tgz", ".tar.xz", ".tar", ".zip"):
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


_CONVERTER_APP = None


def _get_converter_app():
    global _CONVERTER_APP
    if _CONVERTER_APP is None:
        try:
            from flask import Flask
        except Exception:
            return None
        app = Flask("chemotion_converter_app")
        profiles_dir = os.path.join(os.path.dirname(__file__), "chemotion_converter_profiles")
        os.makedirs(os.path.join(profiles_dir, "default"), exist_ok=True)
        app.config["PROFILES_DIR"] = profiles_dir
        _CONVERTER_APP = app
    return _CONVERTER_APP


def _tar_dir_to_temp(path: str) -> Optional[tempfile.NamedTemporaryFile]:
    if not os.path.isdir(path):
        return None
    tf = tempfile.NamedTemporaryFile(suffix=".tar.gz")
    with tarfile.open(tf.name, "w:gz") as tar:
        tar.add(path, arcname=os.path.basename(path))
    tf.seek(0)
    return tf


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


def _normalize_name(value: str) -> str:
    return str(value).strip().lower()


def _column_index(columns: List[Dict], candidates: List[str]) -> Optional[int]:
    wanted = {_normalize_name(c) for c in candidates}
    for idx, col in enumerate(columns):
        name = _normalize_name(col.get("name", ""))
        if name in wanted:
            return idx
    return None


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


def _normalize_tic_label(value) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip().lower()
    if not text:
        return None
    if text in ("plus", "pos", "positive", "+"):
        return "positive"
    if text in ("minus", "neg", "negative", "-"):
        return "negative"
    if text in ("neutral", "neut", "none"):
        return "neutral"
    return None


def _tic_label_from_path(path: Optional[str], tic_hint=None) -> str:
    label = _normalize_tic_label(tic_hint)
    if label:
        return label
    if path:
        name = os.path.basename(path).lower()
        if any(token in name for token in ("plus", "pos", "positive", "positiv")):
            return "positive"
        if any(token in name for token in ("minus", "neg", "negative", "negativ", "minus")):
            return "negative"
        if "neutral" in name or "neut" in name:
            return "neutral"
    return "neutral"


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


def _extract_ms_page(content: str, target_page=None) -> Tuple[List[float], List[float], Optional[float], Optional[str]]:
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


def _extract_xy(table: Dict, x_candidates: List[str], y_candidates: List[str]) -> Tuple[List[float], List[float]]:
    columns = table.get("columns", [])
    rows = table.get("rows", [])
    x_idx = _column_index(columns, x_candidates)
    y_idx = _column_index(columns, y_candidates)
    if x_idx is None or y_idx is None:
        return [], []

    xs: List[float] = []
    ys: List[float] = []
    for row in rows:
        if x_idx >= len(row) or y_idx >= len(row):
            continue
        x_val = _as_float(row[x_idx])
        y_val = _as_float(row[y_idx])
        if x_val is None or y_val is None:
            continue
        xs.append(x_val)
        ys.append(y_val)
    return xs, ys


def _table_metadata(table: Dict) -> Dict[str, str]:
    metadata = table.get("metadata", {}) or {}
    return {str(k).lower(): str(v) for k, v in metadata.items()}


def _resolve_mode(meta: Dict[str, str]) -> Optional[str]:
    mode = meta.get("mode", "").lower()
    if "positiv" in mode:
        return "pos"
    if "negativ" in mode:
        return "neg"
    return None


def _build_ms_rows(table: Dict) -> Tuple[List[Tuple[float, float, float]], Optional[str]]:
    meta = _table_metadata(table)
    mode = _resolve_mode(meta)
    rt_val = meta.get("t") or meta.get("rt") or meta.get("retention_time")
    rt = _as_float(rt_val)
    if rt is None:
        return [], mode

    columns = table.get("columns", [])
    rows = table.get("rows", [])
    mz_idx = _column_index(columns, ["mz", "m/z"])
    int_idx = _column_index(columns, ["intensities", "intensity"])
    if mz_idx is None or int_idx is None:
        return [], mode

    out: List[Tuple[float, float, float]] = []
    for row in rows:
        if mz_idx >= len(row) or int_idx >= len(row):
            continue
        mz = _as_float(row[mz_idx])
        inten = _as_float(row[int_idx])
        if mz is None or inten is None:
            continue
        out.append((rt, mz, inten))
    return out, mode


def lcms_frames_from_converter_app(
    source_path: str,
) -> Optional[Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Optional[str]]]:
    try:
        from converter_app.readers.lcms_reader import LcmsReader  # type: ignore
        from converter_app.readers.helper.reader import Readers  # type: ignore
        from converter_app.models import File as ConverterFile  # type: ignore
    except Exception:
        return None

    if not source_path or not os.path.exists(source_path):
        return None

    class _FakeUpload:
        def __init__(self, path: str):
            self.filename = path
            self.content_type = "application/octet-stream"
            self._fp = open(path, "rb")

        def read(self, *args, **kwargs):
            return self._fp.read(*args, **kwargs)

        def seek(self, *args, **kwargs):
            return self._fp.seek(*args, **kwargs)

        def save(self, dst):
            with open(dst, "wb") as out:
                self._fp.seek(0)
                out.write(self._fp.read())

        def close(self):
            try:
                self._fp.close()
            except Exception:
                pass

    tables = None
    if os.path.isdir(source_path):
        cdf_path = _find_first_cdf(source_path)
        if not cdf_path:
            return None
        file_main = _FakeUpload(cdf_path)
        file_aux = _FakeUpload(cdf_path)
        try:
            main = ConverterFile(file_main)
            aux = ConverterFile(file_aux)
            reader = LcmsReader(main, aux)
            tables = reader.get_tables()
        finally:
            file_main.close()
            file_aux.close()
    elif _is_tar_path(source_path):
        file_main = _FakeUpload(source_path)
        try:
            main = ConverterFile(file_main)
            reader = Readers.instance().match_reader(main)
            if reader is None:
                return None
            tables = reader.get_tables()
        finally:
            file_main.close()
    else:
        return None

    if not tables:
        return None

    uv_rows: List[Dict[str, float]] = []
    ms_neg_rows: List[Tuple[float, float, float]] = []
    ms_pos_rows: List[Tuple[float, float, float]] = []
    ms_generic_rows: List[Tuple[float, float, float]] = []
    has_explicit_mode = False
    has_unknown_mode = False

    for table in tables:
        meta = _table_metadata(table)
        reader_type = meta.get("internal_reader_type", "").lower()

        if reader_type == "lc - uv/vis" or ("wave" in meta or "allwaves" in meta):
            wavelength = meta.get("wavelength") or meta.get("wave") or meta.get("wl")
            if wavelength is None:
                continue
            xs, ys = _extract_xy(
                table,
                ["retentiontime", "retention_time", "time"],
                ["detectorsignal", "detector_signal", "signal", "wavelength"],
            )
            for x, y in zip(xs, ys):
                uv_rows.append(
                    {
                        "RetentionTime": x,
                        "DetectorSignal": y,
                        "wavelength": wavelength,
                    }
                )
            continue

        if "ms spectrum" in reader_type or "mass spectrum" in reader_type:
            rows, mode = _build_ms_rows(table)
            if not rows:
                continue
            if mode == "neg":
                ms_neg_rows.extend(rows)
                has_explicit_mode = True
            elif mode == "pos":
                ms_pos_rows.extend(rows)
                has_explicit_mode = True
            else:
                ms_generic_rows.extend(rows)
                has_unknown_mode = True
            continue

    if not uv_rows and not ms_neg_rows and not ms_pos_rows and not ms_generic_rows:
        return None

    polarity_hint: Optional[str] = None
    if has_unknown_mode and not has_explicit_mode:
        polarity_hint = "generic"

    if ms_generic_rows:
        ms_pos_rows.extend(ms_generic_rows)
    lc_df = pd.DataFrame(uv_rows, columns=["RetentionTime", "DetectorSignal", "wavelength"])
    ms_neg_df = pd.DataFrame(ms_neg_rows, columns=["time", "mz", "intensities"])
    ms_pos_df = pd.DataFrame(ms_pos_rows, columns=["time", "mz", "intensities"])

    return lc_df, ms_neg_df, ms_pos_df, polarity_hint


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

    polarity: Optional[str] = None
    if "NEG" in combined_upper or "MINUS" in combined_upper:
        polarity = "minus"
    elif "POS" in combined_upper or "PLUS" in combined_upper:
        polarity = "plus"
    return kind, polarity


def lcms_preview_image_from_jdx_files(
    jdx_files: List,
    params: Optional[Dict] = None,
) -> Optional[tempfile.NamedTemporaryFile]:
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
    tic_hint = normalized_params.get("lcms_tic")
    mz_page = normalized_params.get("lcms_mz_page")

    peak_path = _pick_jdx_path(paths, "uvvis.peak", "peak.jdx")
    tic_path = _pick_jdx_path(paths, "_tic", "tic")
    mz_path = _pick_jdx_path(paths, "_mz", "mz")

    if tic_hint:
        hint = str(tic_hint).lower()
        for path in paths:
            name = os.path.basename(path).lower()
            if hint in name:
                tic_path = path
                break
        if tic_path is None:
            if hint in ("plus", "pos", "positive", "+"):
                for path in paths:
                    name = os.path.basename(path).lower()
                    if "tic" in name and ("plus" in name or "pos" in name):
                        tic_path = path
                        break
            elif hint in ("minus", "neg", "negative", "-"):
                for path in paths:
                    name = os.path.basename(path).lower()
                    if "tic" in name and ("minus" in name or "neg" in name):
                        tic_path = path
                        break

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

    if not tic_path or not mz_path:
        for path in paths:
            if path == peak_path:
                continue
            try:
                with open(path, "r", encoding="utf-8", errors="ignore") as handle:
                    content = handle.read()
                kind = _classify_lcms_content(content)
                if kind == "tic" and not tic_path:
                    tic_path = path
                elif kind == "mz" and not mz_path:
                    mz_path = path
            except Exception:
                continue
    if not tic_path:
        for path in paths:
            if path == peak_path or path == mz_path:
                continue
            try:
                with open(path, "r", encoding="utf-8", errors="ignore") as handle:
                    content = handle.read()
                if "CHEMSPECTRA UVVIS PEAK TABLE" in content:
                    continue
                xs, ys = _extract_xy_from_jdx_content(content)
                if xs and ys:
                    tic_path = path
                    break
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

    tic_data = None
    if tic_path:
        try:
            with open(tic_path, "r", encoding="utf-8", errors="ignore") as handle:
                tic_content = handle.read()
            tic_data = _extract_xy_from_jdx_content(tic_content)
        except Exception:
            tic_data = None

    ms_data = None
    ms_threshold = None
    ms_page_label = None
    if mz_path:
        try:
            with open(mz_path, "r", encoding="utf-8", errors="ignore") as handle:
                mz_content = handle.read()
            ms_data = _extract_ms_page(mz_content, mz_page)
            if ms_data:
                ms_threshold = ms_data[2]
                ms_page_label = _format_ms_page_label(ms_data[3]) or _format_ms_page_label(mz_page)
        except Exception:
            ms_data = None

    has_uvvis = uvvis_data is not None
    has_tic = tic_data is not None and tic_data[0]
    has_ms = ms_data is not None and ms_data[0]
    if not (has_uvvis or has_tic or has_ms):
        return None

    try:
        import matplotlib.pyplot as plt  # type: ignore
        import matplotlib.path as mpath  # type: ignore
        import numpy as np  # type: ignore
        from chem_spectra.lib.shared.calc import get_curve_endpoint, cal_slope
    except Exception:
        return None

    plt.rcParams["figure.figsize"] = [16, 12]
    plt.rcParams["figure.dpi"] = 200
    plt.rcParams["font.size"] = 14

    fig, axes = plt.subplots(3, 1)
    uvvis_ax, tic_ax, ms_ax = axes

    if has_uvvis:
        xs, ys, edit_peaks, integrations, _wl = uvvis_data
        uvvis_ax.plot(xs, ys)

        if edit_peaks:
            path_data = [
                (mpath.Path.MOVETO, (0, 5)),
                (mpath.Path.LINETO, (0, 20)),
            ]
            codes, verts = zip(*path_data)
            marker = mpath.Path(verts, codes)
            x_peaks = [p["x"] for p in edit_peaks]
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
                x_left, x_right = itg["xL"], itg["xU"]
                i_left, i_right = get_curve_endpoint(xs, ys, x_left, x_right)
                cxs = xs[i_left:i_right]
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

        uvvis_ax.set_xlabel("X (Retention Time)", fontsize=18)
        uvvis_ax.set_ylabel("Y (Detector Signal)", fontsize=18)
        uvvis_ax.grid(False)
        if xs.size:
            x_min, x_max = float(np.min(xs)), float(np.max(xs))
        else:
            x_min, x_max = 0.0, 1.0
        if ys.size:
            y_min, y_max = float(np.min(ys)), float(np.max(ys))
        else:
            y_min, y_max = 0.0, 1.0
        h = max(y_max - y_min, 1.0)
        uvvis_ax.set_xlim(x_min, x_max)
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
    else:
        uvvis_ax.text(0.5, 0.5, "UVVIS unavailable", ha="center", va="center", transform=uvvis_ax.transAxes)
        uvvis_ax.set_axis_off()

    if has_tic:
        tic_xs, tic_ys = tic_data
        tic_ax.plot(tic_xs, tic_ys)
        tic_ax.set_xlabel("X (Retention Time)", fontsize=18)
        tic_ax.set_ylabel("Y (TIC)", fontsize=18)
        tic_ax.grid(False)
        if tic_xs and tic_ys:
            x_min, x_max = min(tic_xs), max(tic_xs)
            y_min, y_max = min(tic_ys), max(tic_ys)
            h = max(y_max - y_min, 1.0)
            tic_ax.set_xlim(x_min, x_max)
            tic_ax.set_ylim(y_min - h * 0.05, y_max + h * 0.15)
        tic_label = _tic_label_from_path(tic_path, tic_hint)
        if tic_label:
            tic_ax.text(
                0.98,
                0.98,
                tic_label,
                ha="right",
                va="top",
                transform=tic_ax.transAxes,
            )
    else:
        tic_ax.text(0.5, 0.5, "TIC unavailable", ha="center", va="center", transform=tic_ax.transAxes)
        tic_ax.set_axis_off()

    if has_ms:
        ms_xs, ms_ys, ms_threshold, _page = ms_data
        threshold = ms_threshold if ms_threshold is not None else 0.05
        max_y = max(ms_ys) if ms_ys else 0.0
        cut = max_y * threshold
        blues_x, blues_y, greys_x, greys_y = [], [], [], []
        for x_val, y_val in zip(ms_xs, ms_ys):
            if y_val >= cut:
                blues_x.append(x_val)
                blues_y.append(y_val)
            else:
                greys_x.append(x_val)
                greys_y.append(y_val)
        ms_ax.bar(greys_x, greys_y, width=0, edgecolor="#dddddd")
        ms_ax.bar(blues_x, blues_y, width=0, edgecolor="#1f77b4")
        ms_ax.set_xlabel("X (m/z)", fontsize=18)
        ms_ax.set_ylabel("Y (Relative Abundance)", fontsize=18)
        ms_ax.grid(False)
        if ms_page_label:
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
    tf_img = tempfile.NamedTemporaryFile(suffix="_lcms_preview.png")
    plt.savefig(tf_img, format="png")
    tf_img.seek(0)
    plt.clf()
    plt.cla()
    plt.close(fig)
    return tf_img


def lcms_jcamp_files_from_converter_app(
    source_path: str,
    title: str,
    lc_df: Optional[pd.DataFrame] = None,
    params: Optional[Dict] = None,
) -> Optional[List[tempfile.NamedTemporaryFile]]:
    try:
        from werkzeug.datastructures import FileStorage
        from converter_app.converters import Converter  # type: ignore
        from converter_app.models import File as ConverterFile  # type: ignore
        from converter_app.readers import READERS  # type: ignore
        from converter_app.writers.jcamp import JcampWriter  # type: ignore
    except Exception:
        return None

    if not source_path or not os.path.exists(source_path):
        return None

    app = _get_converter_app()
    if app is None:
        return None

    tar_temp = None
    src_path = source_path
    if os.path.isdir(source_path):
        tar_temp = _tar_dir_to_temp(source_path)
        if tar_temp is None:
            return None
        src_path = tar_temp.name

    file_handle = None
    try:
        file_handle = open(src_path, "rb")
        fs = FileStorage(
            stream=file_handle,
            filename=src_path,
            content_type="application/octet-stream",
        )
        converter_file = ConverterFile(fs)

        with app.app_context():
            reader = READERS.match_reader(converter_file)
            if not reader:
                return None
            reader.process()
            for table in reader.tables or []:
                header = table.get("header") or []
                if not header:
                    metadata = table.get("metadata", {}) or {}
                    reader_type = metadata.get("internal_reader_type")
                    if reader_type:
                        table["header"] = [str(reader_type)]
                        continue
                    meta_keys = {str(k).lower() for k in metadata.keys()}
                    if {"wave", "wavelength", "allwaves"} & meta_keys:
                        table["header"] = ["lc - uv/vis"]

            converter = Converter.match_profile("default", reader.as_dict)
            if not converter:
                converter = None
            elif converter:
                converter.process()

        jdx_entries: List[Tuple[Dict, bytes]] = []

        if converter is not None:
            jc = JcampWriter(converter)
            for tables in jc.process_ntuples_tables():
                header = tables[0].get("header", {})
                content = jc.write()
                if isinstance(content, str):
                    content = content.encode("utf-8")
                jdx_entries.append((header, content))

            for table in [t for t in converter.tables if t.get("header", {}).get("DATA CLASS") != "NTUPLES"]:
                jc = JcampWriter(converter)
                jc.process_table(table)
                content = jc.write()
                if isinstance(content, str):
                    content = content.encode("utf-8")
                jdx_entries.append((table.get("header", {}), content))

        base = _strip_archive_suffix(title)

        grouped = {"uvvis": [], "tic": [], "mz": []}
        for header, content in jdx_entries:
            kind, polarity = _classify_lcms_header(header)
            if not kind:
                continue
            grouped[kind].append({"polarity": polarity, "content": content})

        output_files: List[tempfile.NamedTemporaryFile] = []

        has_primary = False
        if grouped["uvvis"]:
            output_files.append(_write_named_file(grouped["uvvis"][0]["content"], f"{base}_uvvis.jdx"))
            has_primary = True

        tic_entries = grouped["tic"]
        if len(tic_entries) == 1:
            output_files.append(_write_named_file(tic_entries[0]["content"], f"{base}_tic.jdx"))
            has_primary = True
        elif len(tic_entries) > 1:
            for entry in tic_entries:
                suffix = "_tic.jdx"
                if entry["polarity"] == "plus":
                    suffix = "_tic_plus.jdx"
                elif entry["polarity"] == "minus":
                    suffix = "_tic_minus.jdx"
                output_files.append(_write_named_file(entry["content"], f"{base}{suffix}"))
                has_primary = True

        mz_entries = grouped["mz"]
        if len(mz_entries) == 1:
            output_files.append(_write_named_file(mz_entries[0]["content"], f"{base}_mz.jdx"))
            has_primary = True
        elif len(mz_entries) > 1:
            for entry in mz_entries:
                suffix = "_mz.jdx"
                if entry["polarity"] == "plus":
                    suffix = "_mz_plus.jdx"
                elif entry["polarity"] == "minus":
                    suffix = "_mz_minus.jdx"
                output_files.append(_write_named_file(entry["content"], f"{base}{suffix}"))
                has_primary = True

        if lc_df is None:
            frames = lcms_frames_from_converter_app(source_path)
            lc_df = frames[0] if frames else None

        if lc_df is not None and not lc_df.empty:
            normalized_params = _normalize_params(params)
            uvvis_peak_file = lcms_uvvis_peak_jcamp_from_df(lc_df, title, normalized_params)
            if uvvis_peak_file is not None:
                uvvis_peak_file.seek(0)
                peak_content = uvvis_peak_file.read()
                output_files.append(_write_named_file(peak_content, f"{base}_uvvis.peak.jdx"))

        if output_files and has_primary:
            return output_files
        return None
    finally:
        if file_handle is not None:
            try:
                file_handle.close()
            except Exception:
                pass
        if tar_temp is not None:
            try:
                tar_temp.close()
            except Exception:
                pass


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


def lcms_uvvis_image_from_peak_jdx(peak_jdx_path: str) -> Optional[tempfile.NamedTemporaryFile]:
    """Generate LCMS UVVIS image directly from uvvis.peak.jdx file with peaks and integrations."""
    if not peak_jdx_path or not os.path.exists(peak_jdx_path):
        return None
    
    try:
        import matplotlib.pyplot as plt  # type: ignore
        import matplotlib.path as mpath  # type: ignore
        import numpy as np  # type: ignore
        from chem_spectra.lib.shared.calc import calc_ks, get_curve_endpoint, cal_slope
    except Exception as e:
        return None
    
    try:
        with open(peak_jdx_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        
        # Check if it's a uvvis peak file
        if 'CHEMSPECTRA UVVIS PEAK TABLE' not in content:
            return None
        
        # Collect data by wavelength
        data_by_wavelength = {}
        peaks_by_wavelength = {}
        integrations_by_wavelength = {}
        
        # Split by PAGE blocks in UVVIS PEAK TABLE sections
        sections = content.split('$$ === CHEMSPECTRA UVVIS PEAK TABLE ===')
        for section in sections[1:]:  # Skip first part before first section
            # Extract PAGE value (wavelength)
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
            
            # Extract DATA TABLE section
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
                    # Parse XY data: "x, y;" or "x, y"
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
            
            if xs and ys:
                data_by_wavelength[wavelength] = (xs, ys)
        
        # Extract peaks from EDIT_PEAK sections
        edit_sections = content.split('$$ === CHEMSPECTRA PEAK TABLE EDIT ===')
        for section in edit_sections[1:]:
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
            
            peaks = []
            in_peak_section = False
            for line in section.split('\n'):
                if '##PEAKTABLE=' in line.upper():
                    in_peak_section = True
                    continue
                if '##END=' in line or '$$' in line:
                    break
                if in_peak_section and line.strip() and not line.strip().startswith('##'):
                    parts = line.strip().split(',')
                    if len(parts) >= 2:
                        try:
                            x = _as_float(parts[0].strip())
                            y = _as_float(parts[1].strip())
                            if x is not None and y is not None:
                                peaks.append({'x': x, 'y': y})
                        except (ValueError, IndexError):
                            continue
            
            if peaks:
                peaks_by_wavelength[wavelength] = peaks
        
        # Extract integrations from INTEGRALS sections
        intg_sections = content.split('$$ === CHEMSPECTRA INTEGRALS AND MULTIPLETS ===')
        for section in intg_sections[1:]:
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
            
            integrations = []
            in_intg_section = False
            for line in section.split('\n'):
                if '##$OBSERVEDINTEGRALS=' in line.upper():
                    in_intg_section = True
                    continue
                if '##END=' in line or '$$' in line:
                    break
                if in_intg_section and line.strip() and not line.strip().startswith('##'):
                    # Parse integration: "(xL, xU, area, absoluteArea)"
                    cleaned = line.strip().replace('(', '').replace(')', '')
                    parts = cleaned.split(',')
                    if len(parts) >= 3:
                        try:
                            xL = _as_float(parts[0].strip())
                            xU = _as_float(parts[1].strip())
                            area = _as_float(parts[2].strip())
                            if xL is not None and xU is not None and area is not None:
                                integrations.append({'xL': xL, 'xU': xU, 'area': area})
                        except (ValueError, IndexError):
                            continue
            
            if integrations:
                integrations_by_wavelength[wavelength] = integrations
        
        if not data_by_wavelength:
            return None
        
        # Select wavelength with minimum value for image
        def _to_float_or_none(k):
            try:
                return float(k)
            except Exception:
                return None
        
        numeric_pairs = [(wl, _to_float_or_none(wl)) for wl in data_by_wavelength.keys()]
        numeric_only = [p for p in numeric_pairs if p[1] is not None]
        wl_key = min(numeric_only, key=lambda p: p[1])[0] if numeric_only else list(data_by_wavelength.keys())[0]
        
        xs, ys = data_by_wavelength[wl_key]
        xs = np.asarray(xs, float)
        ys = np.asarray(ys, float)
        
        # Get peaks and integrations for selected wavelength
        edit_peaks = peaks_by_wavelength.get(wl_key, [])
        integrations = integrations_by_wavelength.get(wl_key, [])

        plt.rcParams["figure.figsize"] = [16, 9]
        plt.rcParams["figure.dpi"] = 200
        plt.rcParams["font.size"] = 14

        # Plot main curve
        plt.plot(xs, ys)
        
        # Plot red peaks
        if edit_peaks:
            path_data = [
                (mpath.Path.MOVETO, (0, 5)),
                (mpath.Path.LINETO, (0, 20)),
            ]
            codes, verts = zip(*path_data)
            marker = mpath.Path(verts, codes)
            
            x_peaks = [p['x'] for p in edit_peaks]
            y_peaks = [p['y'] for p in edit_peaks]
            
            plt.plot(
                x_peaks,
                y_peaks,
                'r',
                ls='',
                marker=marker,
                markersize=50,
            )
        
        # Plot integrations (fill area under curve like HPLC UVVIS)
        if integrations:
            y_max = float(np.max(ys))
            y_min = float(np.min(ys))
            h = max(y_max - y_min, 1.0)
            
            for idx, itg in enumerate(integrations):
                xL, xU = itg['xL'], itg['xU']
                # Find curve endpoints
                iL, iU = get_curve_endpoint(xs, ys, xL, xU)
                cxs = xs[iL:iU]
                cys = ys[iL:iU]
                
                if len(cxs) > 0 and len(cys) > 0:
                    # Fill area under curve like HPLC UVVIS
                    slope = cal_slope(cxs[0], cys[0], cxs[len(cxs)-1], cys[len(cys)-1])
                    last_y = cys[0]
                    last_x = cxs[0]
                    aucys = [last_y]
                    for i in range(1, len(cys)):
                        curr_x = cxs[i]
                        curr_y = slope*(curr_x-last_x) + last_y
                        aucys.append(curr_y)
                        last_x = curr_x
                        last_y = curr_y
                    plt.fill_between(cxs, y1=cys, y2=aucys, alpha=0.2, color='#FF0000')
        
        plt.xlabel("X (Retention Time)", fontsize=18)
        plt.ylabel("Y (Detector Signal)", fontsize=18)
        plt.grid(False)

        if xs.size:
            x_min, x_max = float(np.min(xs)), float(np.max(xs))
        else:
            x_min, x_max = 0.0, 1.0
        if ys.size:
            y_min, y_max = float(np.min(ys)), float(np.max(ys))
        else:
            y_min, y_max = 0.0, 1.0
        h = max(y_max - y_min, 1.0)
        plt.xlim(x_min, x_max)
        plt.ylim(y_min - h * 0.05, y_max + h * 0.15)

        tf_img = tempfile.NamedTemporaryFile(suffix="_lcms_uvvis.peak.png")
        plt.savefig(tf_img, format="png")
        tf_img.seek(0)
        plt.clf()
        plt.cla()
        plt.close()
        return tf_img
    except Exception as e:  # noqa: E722
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
            f"##PAGE={wl}\n",
            f"##NPOINTS={len(xs)}\n",
            "##DATA TABLE= (XY..XY), PEAKS\n",
        ]
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

