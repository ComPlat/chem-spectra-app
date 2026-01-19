import logging
import os
from typing import Dict, List, Optional, Tuple

import pandas as pd

logger = logging.getLogger(__name__)


def _find_first_cdf(root_dir: str) -> Optional[str]:
    for dirpath, _, filenames in os.walk(root_dir):
        for fn in filenames:
            if fn.lower().endswith(".cdf"):
                return os.path.join(dirpath, fn)
    return None


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


def lcms_frames_from_converter_app(root_dir: str) -> Optional[Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]]:
    """
    Lit un dossier LCMS via le lecteur chemotion-converter-app (lcms_reader).
    Retourne (lc_df, ms_minus_df, ms_plus_df) si possible, sinon None.
    """
    try:
        from converter_app.readers.lcms_reader import LcmsReader  # type: ignore
        from converter_app.models import File as ConverterFile  # type: ignore
    except Exception:
        return None

    cdf_path = _find_first_cdf(root_dir)
    if not cdf_path:
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

    file_main = _FakeUpload(cdf_path)
    file_aux = _FakeUpload(cdf_path)
    try:
        main = ConverterFile(file_main)
        aux = ConverterFile(file_aux)
        reader = LcmsReader(main, aux)
        tables = reader.prepare_tables()
    finally:
        file_main.close()
        file_aux.close()

    if not tables:
        return None

    uv_rows: List[Dict[str, float]] = []
    ms_neg_rows: List[Tuple[float, float, float]] = []
    ms_pos_rows: List[Tuple[float, float, float]] = []
    tic_neg: List[Tuple[float, float]] = []
    tic_pos: List[Tuple[float, float]] = []

    for table in tables:
        meta = _table_metadata(table)
        reader_type = meta.get("internal_reader_type", "").lower()

        if reader_type == "lc - uv/vis":
            wavelength = meta.get("wavelength") or meta.get("wl")
            if wavelength is None:
                continue
            xs, ys = _extract_xy(
                table,
                ["retentiontime", "retention_time", "time"],
                ["detectorsignal", "detector_signal", "signal"],
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

        if reader_type == "ms spectrum":
            rows, mode = _build_ms_rows(table)
            if not rows:
                continue
            if mode == "neg":
                ms_neg_rows.extend(rows)
            else:
                ms_pos_rows.extend(rows)
            continue

        if reader_type == "ms chromatogramm":
            xs, ys = _extract_xy(
                table,
                ["time", "retention_time"],
                ["tic", "intensity", "intensities"],
            )
            mode = _resolve_mode(meta)
            if mode == "neg":
                tic_neg.extend(zip(xs, ys))
            else:
                tic_pos.extend(zip(xs, ys))

    if not uv_rows and not ms_neg_rows and not ms_pos_rows:
        return None

    lc_df = pd.DataFrame(uv_rows, columns=["RetentionTime", "DetectorSignal", "wavelength"])
    ms_neg_df = pd.DataFrame(ms_neg_rows, columns=["time", "mz", "intensities"])
    ms_pos_df = pd.DataFrame(ms_pos_rows, columns=["time", "mz", "intensities"])

    # Si les TIC ne sont pas fournis, ils seront recalculés plus loin.
    if not tic_neg and not ms_neg_df.empty:
        pass
    if not tic_pos and not ms_pos_df.empty:
        pass

    return lc_df, ms_neg_df, ms_pos_df

