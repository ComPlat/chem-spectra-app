import tempfile
import json
import logging
from typing import List, Optional, Dict

from chem_spectra.lib.external.chemotion_converter_lcms import (
    lcms_preview_image_from_jdx_files,
    lcms_df_from_peak_jdx,
    lcms_df_from_uvvis_jdx,
    lcms_uvvis_peak_jcamp_from_df,
    lcms_uvvis_units_from_jdx,
)

logger = logging.getLogger(__name__)

UVVIS_PEAK_MARKER = "CHEMSPECTRA UVVIS PEAK TABLE"

_UVVIS_DATA_TYPE_TOKENS = (
    "HPLC UV-VIS",
    "HPLC UV/VIS",
    "UV-VIS",
    "UV/VIS",
    "ULTRAVIOLET",
)

_MARKER_SCAN_CHUNK = 65536


def has_uvvis_peak_marker(path: Optional[str]) -> bool:
    if not path:
        return False
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as handle:
            while True:
                chunk = handle.read(_MARKER_SCAN_CHUNK)
                if not chunk:
                    return False
                if UVVIS_PEAK_MARKER in chunk:
                    return True
    except OSError:
        return False


class LCMSConverterAppComposer:
    def __init__(
        self,
        jcamp_files: List[tempfile.NamedTemporaryFile],
        image: Optional[tempfile.NamedTemporaryFile] = None,
        params: Optional[Dict] = None,
    ):
        self.data = jcamp_files
        self._image = image
        self.params = params
        self._peaks_applied = False
        self._ensure_uvvis_peak_file()

    @staticmethod
    def _has_lcms_edits(params: Optional[Dict]) -> bool:
        if not params:
            return False
        peaks = params.get("peaks_str")
        if peaks not in (None, "", [], {}):
            return True

        integration = params.get("integration")
        if integration in (None, "", "{}", "[]", {}, []):
            return False
        if isinstance(integration, str):
            try:
                integration = json.loads(integration)
            except Exception:
                return True
        if isinstance(integration, dict):
            stack = integration.get("stack") if "stack" in integration else None
            curves = integration.get("curves") or integration.get("byCurve") or integration.get("by_curve")
            if stack:
                return True
            if isinstance(curves, list) and len(curves) > 0:
                return True
            # Wavelength keyed integrations: {"230": [[xL,xU,area,...]]}
            for value in integration.values():
                if isinstance(value, list) and len(value) > 0:
                    return True
            return False
        return bool(integration)

    @staticmethod
    def _requested_lcms_mz_page(params: Optional[Dict]) -> str:
        if not params:
            return ""
        raw_value = params.get("lcms_mz_page")
        if raw_value is None:
            return ""
        return " ".join(str(raw_value).splitlines()).strip()

    @staticmethod
    def _read_lcms_mz_page_from_file(jdx_file: Optional[tempfile.NamedTemporaryFile]) -> str:
        if not jdx_file:
            return ""
        try:
            with open(jdx_file.name, "r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    if line.startswith("##$CSLCMSMZPAGE="):
                        return line.split("=", 1)[1].strip()
        except Exception:
            return ""
        return ""

    def _should_refresh_jcamp(self) -> bool:
        if self._has_lcms_edits(self.params):
            return True

        requested_mz_page = self._requested_lcms_mz_page(self.params)
        if not requested_mz_page:
            return False

        jdx_file = self._find_peak_or_edit_file()
        if not jdx_file and self.data:
            jdx_file = self.data[0]

        current_mz_page = self._read_lcms_mz_page_from_file(jdx_file)
        if current_mz_page != requested_mz_page:
            return True
        return False

    def tf_img(self):
        if self._should_refresh_jcamp() and not self._peaks_applied:
            self.tf_jcamp()
        if self._image is not None:
            return self._image
        if self.data:
            preview = lcms_preview_image_from_jdx_files(self.data, self.params)
            if preview:
                self._image = preview
        return self._image

    def _find_peak_or_edit_file(self) -> Optional[tempfile.NamedTemporaryFile]:
        for jdx_file in self.data:
            if jdx_file.name.lower().endswith(('peak.jdx', 'edit.jdx')):
                return jdx_file
        return None

    @staticmethod
    def _is_uvvis_source(jdx_file) -> bool:
        path = getattr(jdx_file, 'name', None) if jdx_file else None
        if not path:
            return False
        try:
            with open(path, 'r', encoding='utf-8', errors='ignore') as handle:
                for _ in range(200):
                    raw = handle.readline()
                    if not raw:
                        break
                    line = raw.strip()
                    if not line.startswith('##'):
                        if line:
                            continue
                        break
                    if not line.upper().startswith('##DATA TYPE'):
                        continue
                    value = line.split('=', 1)[1].strip().upper() if '=' in line else ''
                    return any(token in value for token in _UVVIS_DATA_TYPE_TOKENS)
        except OSError:
            return False
        return False

    def _resolve_uvvis_title(self) -> str:
        title = ''
        if isinstance(self.params, dict):
            title = self.params.get('fname') or ''
        if title:
            title = title.replace('.jdx', '').replace('.edit', '').replace('.peak', '')
        return title or 'lc_ms_spectrum'

    def _ensure_uvvis_peak_file(self) -> None:
        if not self.data:
            return

        for jdx_file in self.data:
            name = (getattr(jdx_file, 'name', '') or '').lower()
            if name.endswith(('peak.jdx', 'edit.jdx')):
                return
            if has_uvvis_peak_marker(getattr(jdx_file, 'name', None)):
                return

        for idx, jdx_file in enumerate(self.data):
            path = getattr(jdx_file, 'name', None)
            if not path or not self._is_uvvis_source(jdx_file):
                continue
            try:
                lc_df = lcms_df_from_uvvis_jdx(path)
            except Exception as err:
                logger.debug("UVVIS peak prebake: parse failed on %s: %s", path, err)
                continue
            if lc_df is None or lc_df.empty:
                continue

            try:
                x_units, y_units = lcms_uvvis_units_from_jdx(path)
                new_file = lcms_uvvis_peak_jcamp_from_df(
                    lc_df,
                    self._resolve_uvvis_title(),
                    self.params,
                    x_units=x_units,
                    y_units=y_units,
                )
            except Exception as err:
                logger.warning(
                    "UVVIS peak prebake: generator raised on %s: %s", path, err,
                )
                return

            if new_file is None:
                return

            self.data.insert(idx + 1, new_file)
            self._image = None
            return

    def tf_jcamp(self):
        if not self.data:
            return None
        
        if self._should_refresh_jcamp() and not self._peaks_applied:
            jdx_file = self._find_peak_or_edit_file()
            if not jdx_file and self.data:
                jdx_file = self.data[0]
            
            if jdx_file:
                lc_df = lcms_df_from_peak_jdx(jdx_file.name)
                units_path = jdx_file.name
                if lc_df is None or lc_df.empty:
                    lc_df = lcms_df_from_uvvis_jdx(jdx_file.name)
                if (lc_df is None or lc_df.empty) and self.data:
                    for candidate in self.data:
                        if candidate is jdx_file:
                            continue
                        maybe_df = lcms_df_from_uvvis_jdx(candidate.name)
                        if maybe_df is not None and not maybe_df.empty:
                            lc_df = maybe_df
                            units_path = candidate.name
                            break
                if lc_df is not None and not lc_df.empty:
                    title = self.params.get('fname', 'lc_ms_spectrum')
                    if title:
                        title = title.replace('.jdx', '').replace('.edit', '').replace('.peak', '')
                    else:
                        title = 'lc_ms_spectrum'
                    x_units, y_units = lcms_uvvis_units_from_jdx(units_path)
                    if not x_units or not y_units:
                        for candidate in self.data:
                            candidate_path = getattr(candidate, 'name', None)
                            if not candidate_path or candidate_path == units_path:
                                continue
                            if not self._is_uvvis_source(candidate):
                                continue
                            src_x, src_y = lcms_uvvis_units_from_jdx(candidate_path)
                            x_units = x_units or src_x
                            y_units = y_units or src_y
                            if x_units and y_units:
                                break
                    new_file = lcms_uvvis_peak_jcamp_from_df(
                        lc_df,
                        title,
                        self.params,
                        x_units=x_units,
                        y_units=y_units,
                    )
                    if new_file:
                        replaced = False
                        for idx, existing in enumerate(self.data):
                            if existing is jdx_file:
                                self.data[idx] = new_file
                                replaced = True
                                break
                        if not replaced:
                            self.data.insert(0, new_file)
                        preview = lcms_preview_image_from_jdx_files(self.data, self.params)
                        if preview:
                            self._image = preview
                        self._peaks_applied = True
                        return new_file
            self._peaks_applied = True
        
        return self.data[0] if self.data else None

    def tf_csv(self):
        return None
