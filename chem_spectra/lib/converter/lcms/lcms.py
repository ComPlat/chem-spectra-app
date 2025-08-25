from __future__ import annotations

import json
import logging
import re
from json import JSONDecodeError
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np

from chem_spectra.lib.converter.lcms.base import LCMSBaseConverter
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter

logger = logging.getLogger(__name__)

SECTION_TIC_POS = "TIC POSITIVE SPECTRUM"
SECTION_TIC_NEG = "TIC NEGATIVE SPECTRUM"
SECTION_UVVIS = "UVVIS SPECTRUM"
SECTION_UVVIS_HPLC = "HPLC UV/VIS SPECTRUM"
SECTION_MZ_POS = "MZ POSITIVE SPECTRUM"
SECTION_MZ_NEG = "MZ NEGATIVE SPECTRUM"

DATA_TAGS_START = ("##DATA TABLE=", "##DATATABLE=", "##XYDATA=", "##XYPOINTS=")
DATA_TAGS_STOP = ("##PEAKTABLE=", "##NPOINTS=")
OBSERVED_INTEGRALS_TAGS = ("##$OBSERVEDINTEGRALS=", "##OBSERVEDINTEGRALS=")

NUMERIC_LINE_RE = re.compile(r"^\s*[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?")

SPLIT_RE = re.compile(r"[ ,;\t]+")


def _clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def _coerce_threshold(params: Optional[dict], default_ratio: float = 0.5) -> float:
    if not params:
        return _clamp01(default_ratio)
    raw = params.get("thres", None)
    if raw is None:
        return _clamp01(default_ratio)
    try:
        val = float(raw)
    except (TypeError, ValueError):
        logger.debug("Invalid 'thres' value %r; fallback to default %s", raw, default_ratio)
        return _clamp01(default_ratio)
    if val > 1.0:
        val = val / 100.0
    return _clamp01(val)


def _parse_float(s: str) -> Optional[float]:
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def _first_number_in_string(s) -> Optional[float]:
    if s is None:
        return None
    try:
        text = str(s)
    except Exception:
        return None
    if not text:
        return None
    m = re.search(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", text)
    if not m:
        return None
    return _parse_float(m.group(0))


def _normalize_numeric_key(key) -> str:
    f = _first_number_in_string(key)
    return str(f) if f is not None else ("" if key is None else str(key))


def _split_numbers(line: str) -> List[float]:
    parts = SPLIT_RE.split(line.strip())
    out: List[float] = []
    for p in parts:
        if p == "":
            continue
        v = _parse_float(p)
        if v is not None:
            out.append(v)
    return out


def _read_text_file(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except Exception:
        try:
            return path.read_text(encoding="latin-1")
        except Exception:
            return path.read_text(encoding="utf-8", errors="replace")

class LCMSConverter:
    def __init__(self, jcamp_or_path, params: Optional[dict] = None, fname: str = "chromatogram.jdx"):
        self.params = params or {}
        self.fname = Path(fname).name
        self.jcamp_content: Optional[str] = None

        self.tic_pos: Dict[str, List[float]] = {"time": [], "intensity": []}
        self.tic_neg: Dict[str, List[float]] = {"time": [], "intensity": []}
        self.uvvis: Dict[str, Dict[str, List]] = {}
        self.mz_pos: Dict[str, Dict[str, List[float]]] = {}
        self.mz_neg: Dict[str, Dict[str, List[float]]] = {}
        self._edit_peaks_from_jcamp: Dict[str, List[Dict[str, float]]] = {}

        if isinstance(jcamp_or_path, JcampBaseConverter):
            self.jcamp_content = getattr(jcamp_or_path, "raw_content", None)
            if not self.params:
                self.params = jcamp_or_path.params or {}
            self.fname = Path(getattr(jcamp_or_path, "fname", self.fname)).name

        elif isinstance(jcamp_or_path, (str, Path)):
            p = Path(str(jcamp_or_path))
            if p.is_dir():
                try:
                    jdx_file = next(f for f in p.iterdir() if f.is_file() and f.suffix.lower() in (".jdx", ".dx", ".jcamp"))
                except StopIteration:
                    jdx_file = None
                if jdx_file:
                    try:
                        self.jcamp_content = _read_text_file(jdx_file)
                        self.fname = jdx_file.name
                    except Exception as exc:
                        logger.warning("Failed reading JCAMP %s: %s", jdx_file, exc)
                        self.jcamp_content = None
                else:
                    logger.info("No JCAMP file found in directory %s.", p)

            elif p.is_file() and p.suffix.lower() in (".jdx", ".dx", ".jcamp"):
                try:
                    self.jcamp_content = _read_text_file(p)
                    self.fname = p.name
                except Exception as exc:
                    logger.warning("Failed reading JCAMP %s: %s", p, exc)
                    self.jcamp_content = None

            else:
                s = str(jcamp_or_path)
                if "##JCAMP-DX" in s or "##TITLE=" in s:
                    self.jcamp_content = s

        if self.jcamp_content:
            self._parse_jcamp_single_pass(self.jcamp_content)
        else:
            try:
                base = LCMSBaseConverter(jcamp_or_path, params=self.params, fname=self.fname)  # type: ignore
                if getattr(base, "data", None):
                    tp, tn, uv, mzp, mzn = base.data
                    self.tic_pos = {"time": tp.get("time", []), "intensity": tp.get("Intensity", [])}
                    self.tic_neg = {"time": tn.get("time", []), "intensity": tn.get("Intensity", [])}
                    self.uvvis = uv or {}
                    self.mz_pos = mzp or {}
                    self.mz_neg = mzn or {}
                logger.info("Initialized LCMSBaseConverter fallback.")
            except Exception as exc:
                logger.info("Skipped LCMSBaseConverter fallback (%s).", exc)

        self._finalize_series()
        if self.jcamp_content:
            try:
                self._parse_edit_peaks_from_jcamp(self.jcamp_content)
            except Exception:
                pass

        self.threshold = _coerce_threshold(self.params, default_ratio=0.5)


        self.edit_peaks = self._process_uvvis_peaks_from_params()
        if not self.edit_peaks and self._edit_peaks_from_jcamp:
            self.edit_peaks = self._edit_peaks_from_jcamp
        self.edit_integrals = self._process_uvvis_integrals_from_params()

        self.boundary = self._make_boundary_from_both(self.tic_pos, self.tic_neg)

        self.data = [
            {"time": self.tic_pos["time"], "Intensity": self.tic_pos["intensity"]},
            {"time": self.tic_neg["time"], "Intensity": self.tic_neg["intensity"]},
            self.uvvis,
            self.mz_pos,
            self.mz_neg,
        ]

    def _parse_jcamp_single_pass(self, content: str) -> None:
        in_section: Optional[str] = None
        read_xy: bool = False

        current_uvvis_wl_key: Optional[str] = None
        current_mz_time_key: Optional[str] = None

        for raw_line in content.splitlines():
            line = raw_line.rstrip("\n\r")
            uline = line.upper().strip()

            if SECTION_TIC_POS in uline:
                in_section, read_xy = "tic_pos", False
                current_uvvis_wl_key = None
                current_mz_time_key = None
                continue
            if SECTION_TIC_NEG in uline:
                in_section, read_xy = "tic_neg", False
                current_uvvis_wl_key = None
                current_mz_time_key = None
                continue
            if (SECTION_UVVIS in uline) or (SECTION_UVVIS_HPLC in uline):
                in_section, read_xy = "uvvis", False
                current_mz_time_key = None
                continue
            if SECTION_MZ_POS in uline:
                in_section, read_xy = "mz_pos", False
                current_uvvis_wl_key = None
                continue
            if SECTION_MZ_NEG in uline:
                in_section, read_xy = "mz_neg", False
                current_uvvis_wl_key = None
                continue

            if in_section is None:
                continue

            if uline.startswith("##PAGE="):
                page_val = line.split("=", 1)[1].strip()
                if in_section == "uvvis":
                    wl = _first_number_in_string(page_val)
                    wl_key = _normalize_numeric_key(page_val if wl is None else str(wl))
                    current_uvvis_wl_key = wl_key
                    self.uvvis.setdefault(wl_key, {"RetentionTime": [], "DetectorSignal": [], "integrals": []})
                    read_xy = False
                    continue
                elif in_section in ("mz_pos", "mz_neg"):
                    rt = _first_number_in_string(page_val)
                    rt_key = _normalize_numeric_key(page_val if rt is None else str(rt))
                    current_mz_time_key = rt_key
                    tgt = self.mz_pos if in_section == "mz_pos" else self.mz_neg
                    tgt.setdefault(rt_key, {"mz": [], "intensities": []})
                    read_xy = False
                continue

            if uline.startswith(DATA_TAGS_START):
                read_xy = True
                continue
            if uline.startswith(DATA_TAGS_STOP) or uline.startswith("$$ === CHEMSPECTRA INTEGRALS"):
                read_xy = False
                continue
            
            if any(uline.startswith(t) for t in OBSERVED_INTEGRALS_TAGS):
                integrals_str = line.split("=", 1)[1].strip()
                tuples = re.findall(r"\((.*?)\)", integrals_str)
                values: List[Tuple[float, ...]] = []
                for tp in tuples:
                    nums = [_parse_float(x.strip()) for x in tp.split(",")]
                    if all(v is not None for v in nums):
                        values.append(tuple(v for v in nums if v is not None))  # type: ignore
                if in_section == "uvvis" and current_uvvis_wl_key:
                    self.uvvis.setdefault(current_uvvis_wl_key, {"RetentionTime": [], "DetectorSignal": [], "integrals": []})
                    self.uvvis[current_uvvis_wl_key]["integrals"].extend(values)
                else:
                    self.uvvis.setdefault("default", {"RetentionTime": [], "DetectorSignal": [], "integrals": []})
                    self.uvvis["default"]["integrals"].extend(values)
                continue

            if uline.startswith("##END="):
                in_section = None
                read_xy = False
                current_uvvis_wl_key = None
                current_mz_time_key = None
                continue

            if read_xy and NUMERIC_LINE_RE.match(line):
                nums = _split_numbers(line)
                if len(nums) < 2:
                            continue
                
                if in_section == "tic_pos":
                    self.tic_pos["time"].append(nums[0])
                    self.tic_pos["intensity"].append(nums[1])
                elif in_section == "tic_neg":
                    self.tic_neg["time"].append(nums[0])
                    self.tic_neg["intensity"].append(nums[1])
                elif in_section == "uvvis":
                    wl_key = current_uvvis_wl_key or "default"
                    self.uvvis.setdefault(wl_key, {"RetentionTime": [], "DetectorSignal": [], "integrals": []})
                    self.uvvis[wl_key]["RetentionTime"].append(nums[0])
                    self.uvvis[wl_key]["DetectorSignal"].append(nums[1])
                elif in_section in ("mz_pos", "mz_neg"):
                    time_key = current_mz_time_key or "0"
                    tgt = self.mz_pos if in_section == "mz_pos" else self.mz_neg
                    tgt.setdefault(time_key, {"mz": [], "intensities": []})
                    tgt[time_key]["mz"].append(nums[0])
                    tgt[time_key]["intensities"].append(nums[1])

    def _finalize_series(self) -> None:
        def _clean_and_sort_xy(x: List[float], y: List[float]) -> Tuple[List[float], List[float]]:
            if not x or not y:
                return [], []
            xa = np.asarray(x, dtype=float)
            ya = np.asarray(y, dtype=float)
            good = np.isfinite(xa) & np.isfinite(ya) 
            xa, ya = xa[good], ya[good]
            if xa.size == 0:
                return [], []
            order = np.argsort(xa)
            return xa[order].tolist(), ya[order].tolist()

        self.tic_pos["time"], self.tic_pos["intensity"] = _clean_and_sort_xy(
            self.tic_pos["time"], self.tic_pos["intensity"]
        )
        self.tic_neg["time"], self.tic_neg["intensity"] = _clean_and_sort_xy(
            self.tic_neg["time"], self.tic_neg["intensity"]
        )

        for wl_key, serie in list(self.uvvis.items()):
            rt = serie.get("RetentionTime", []) or []
            ds = serie.get("DetectorSignal", []) or []
            x, y = _clean_and_sort_xy(rt, ds)
            self.uvvis[wl_key]["RetentionTime"] = x
            self.uvvis[wl_key]["DetectorSignal"] = y

        def _clean_mz_dict(d: Dict[str, Dict[str, List[float]]]) -> None:
            for k, serie in d.items():
                mzs = serie.get("mz", []) or []
                ints = serie.get("intensities", []) or []
                x, y = _clean_and_sort_xy(mzs, ints)
                serie["mz"] = x
                serie["intensities"] = y

        _clean_mz_dict(self.mz_pos)
        _clean_mz_dict(self.mz_neg)

    def _make_boundary_from_both(self, s1: Dict[str, List[float]], s2: Dict[str, List[float]]) -> Dict[str, Dict[str, float]]:
        xs = (s1.get("time") or []) + (s2.get("time") or [])
        ys = (s1.get("intensity") or []) + (s2.get("intensity") or [])
        if not xs or not ys:
            return {"x": {"min": 0.0, "max": 0.0}, "y": {"min": 0.0, "max": 0.0}}
        xa = np.asarray(xs, dtype=float)
        ya = np.asarray(ys, dtype=float)
        good = np.isfinite(xa) & np.isfinite(ya)
        if not np.any(good):
            return {"x": {"min": 0.0, "max": 0.0}, "y": {"min": 0.0, "max": 0.0}}
        xa, ya = xa[good], ya[good]
        return {"x": {"min": float(np.min(xa)), "max": float(np.max(xa))},
                "y": {"min": float(np.min(ya)), "max": float(np.max(ya))}}

    def _process_uvvis_peaks_from_params(self) -> Dict[str, List[Dict[str, float]]]:
        peaks_str = self.params.get("peaks_str")
        if not peaks_str:
            return {}

        try:
            if isinstance(peaks_str, dict):
                return self._normalize_uvvis_key_dict(peaks_str)
            if isinstance(peaks_str, str):
                s = peaks_str.strip()
                if s.startswith("{"):
                    obj = json.loads(s)
                    if isinstance(obj, dict):
                        return self._normalize_uvvis_key_dict(obj)
                    return {}
                pairs = [p for p in s.split("#") if "," in p]
                pts: List[Dict[str, float]] = []
                for p in pairs:
                    parts = p.split(",", 1)
                    if len(parts) == 2:
                        x = _parse_float(parts[0])
                        y = _parse_float(parts[1])
                        if x is not None and y is not None:
                            pts.append({"x": x, "y": y})
                if pts:
                    return {"default": pts}
            return {}
        except JSONDecodeError as exc:
            logger.warning("[UVVIS peaks] JSON parse error: %s", exc)
            return {}
        except Exception as exc:
            logger.warning("[UVVIS peaks] parse error: %s", exc)
            return {}

    def _process_uvvis_integrals_from_params(self) -> Dict[str, List[Tuple[float, ...]]]:
        integration_data = self.params.get("integration")
        if not integration_data:
            return {}
        try:
            if isinstance(integration_data, dict):
                raw = integration_data
            elif isinstance(integration_data, str):
                raw = json.loads(integration_data)
            else:
                return {}

            if isinstance(raw, dict) and "stack" in raw:
                return {}

            if isinstance(raw, dict):
                return self._normalize_uvvis_key_dict(raw)
            return {}
        except JSONDecodeError as exc:
            logger.warning("[UVVIS integrals] JSON parse error: %s", exc)
            return {}
        except Exception as exc:
            logger.warning("[UVVIS integrals] parse error: %s", exc)
            return {}

    @staticmethod
    def _normalize_uvvis_key_dict(d: Dict) -> Dict:
        out = {}
        for k, v in d.items():
            nk = _normalize_numeric_key(str(k))
            out[nk] = v
        return out

    def _parse_edit_peaks_from_jcamp(self, content: str) -> None:
        in_uvvis = False
        collecting = False
        current_wl_key: Optional[str] = None
        for raw_line in content.splitlines():
            line = raw_line.rstrip("\n\r")
            uline = line.upper().strip()

            if 'UVVIS SPECTRUM' in uline or 'HPLC UV/VIS SPECTRUM' in uline:
                in_uvvis = True
                collecting = False
                current_wl_key = None
                continue
            if not in_uvvis:
                continue

            if uline.startswith('##PAGE='):
                page_val = line.split('=', 1)[1].strip()
                wl = _first_number_in_string(page_val)
                wl_key = _normalize_numeric_key(page_val if wl is None else str(wl))
                current_wl_key = wl_key
                collecting = False
                continue

            if uline.startswith('##PEAKTABLE='):
                collecting = True
                # ensure key
                if current_wl_key is None:
                    current_wl_key = 'default'
                continue

            if uline.startswith('##') or uline.startswith('$$') or uline.startswith('##END='):
                collecting = False
                continue

            if collecting and current_wl_key is not None:
                m = re.match(r"^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)[ ,;]+([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)[,;]?\s*;?\s*$", line)
                if m:
                    x = _parse_float(m.group(1))
                    y = _parse_float(m.group(2))
                    if x is not None and y is not None:
                        self._edit_peaks_from_jcamp.setdefault(current_wl_key, []).append({"x": x, "y": y})

    def _parse_jcamp(self):
        if not self.jcamp_content:
            return [
                {"time": [], "Intensity": []},
                {"time": [], "Intensity": []},
                {},
                {},
                {},
            ]
        self._parse_jcamp_single_pass(self.jcamp_content)
        self._finalize_series()
        return [
            {"time": self.tic_pos["time"], "Intensity": self.tic_pos["intensity"]},
            {"time": self.tic_neg["time"], "Intensity": self.tic_neg["intensity"]},
            self.uvvis,
            self.mz_pos,
            self.mz_neg,
        ]

    def _process_uvvis_peaks(self):
        return self._process_uvvis_peaks_from_params()

    def _process_uvvis_integrals(self):
        return self._process_uvvis_integrals_from_params()
