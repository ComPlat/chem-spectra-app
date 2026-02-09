import json
import os
import tarfile
import tempfile
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
    combined = " ".join([v for v in [category, ntuples_id, data_type] if v])
    combined_upper = combined.upper()

    kind: Optional[str] = None
    if "UVVIS" in combined_upper or "UV/VIS" in combined_upper or "UV-VIS" in combined_upper or "UV VIS" in combined_upper:
        kind = "uvvis"
    elif "TIC" in combined_upper:
        kind = "tic"
    elif "MZ" in combined_upper or "M/Z" in combined_upper or "MASS SPECTRUM" in combined_upper:
        kind = "mz"
    if not kind:
        return None, None

    polarity: Optional[str] = None
    if "NEG" in combined_upper or "MINUS" in combined_upper:
        polarity = "minus"
    elif "POS" in combined_upper or "PLUS" in combined_upper:
        polarity = "plus"
    return kind, polarity


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


def lcms_uvvis_image_from_df(lc_df: pd.DataFrame) -> Optional[tempfile.NamedTemporaryFile]:
    if lc_df is None or lc_df.empty:
        return None
    try:
        import matplotlib.pyplot as plt  # type: ignore
        import numpy as np  # type: ignore
    except Exception:
        return None

    wl_key = None
    if "wavelength" in lc_df.columns and not lc_df["wavelength"].isna().all():
        def _to_float_or_none(k):
            try:
                return float(k)
            except Exception:
                return None
        numeric_pairs = [(k, _to_float_or_none(k)) for k in lc_df["wavelength"].unique().tolist()]
        numeric_only = [p for p in numeric_pairs if p[1] is not None]
        wl_key = min(numeric_only, key=lambda p: p[1])[0] if numeric_only else lc_df["wavelength"].iloc[0]

    if wl_key is not None:
        data = lc_df[lc_df["wavelength"] == wl_key]
    else:
        data = lc_df

    xs = np.asarray(data.get("RetentionTime", []), float)
    ys = np.asarray(data.get("DetectorSignal", []), float)

    plt.rcParams["figure.figsize"] = [16, 9]
    plt.rcParams["figure.dpi"] = 200
    plt.rcParams["font.size"] = 14

    plt.plot(xs, ys)
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
                x_left = float(x_left) - float(ref_shift)
                x_right = float(x_right) - float(ref_shift)
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
        integration_meta, integration_stack = _integration_for_page(integration_input, page_idx, wl)
        integration_lines = _integration_lines(integration_meta or integration_input, integration_stack)

        page_block = [
            "\n",
            "$$ === CHEMSPECTRA UVVIS PEAK TABLE ===\n",
            f"##TITLE={title}\n",
            "##JCAMP-DX=5.00\n",
            "##DATA TYPE=LC/MS\n",
            "##DATA CLASS=PEAK TABLE\n",
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

        content.extend([
            "\n",
            "$$ === CHEMSPECTRA PEAK TABLE AUTO ===\n",
            f"##TITLE={title}\n",
            "##JCAMP-DX=5.00\n",
            "##DATA TYPE=LC/MSPEAKTABLE\n",
            "##DATA CLASS=PEAKTABLE\n",
            "##$CSCATEGORY=AUTO_PEAK\n",
            "##$CSTHRESHOLD=0.05\n",
            f"##MAXX={x_max}\n",
            f"##MAXY={y_max}\n",
            f"##MINX={x_min}\n",
            f"##MINY={y_min}\n",
            f"##PAGE={wl}\n",
            f"##NPOINTS={len(auto_peaks)}\n",
            "##PEAKTABLE= (XY..XY)\n",
        ])
        for peak in auto_peaks:
            content.append(f"{peak['x']}, {peak['y']}\n")
        content.extend(["##END=\n"])

    tf = tempfile.NamedTemporaryFile(suffix="_lcms_uvvis.peak.jdx")
    tf.write(bytes("".join(content), "utf-8"))
    tf.seek(0)
    return tf

