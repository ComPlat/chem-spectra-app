import tempfile
import json
from typing import List, Optional, Dict

from chem_spectra.lib.external.chemotion_converter_lcms import (
    lcms_preview_image_from_jdx_files,
    lcms_df_from_peak_jdx,
    lcms_df_from_uvvis_jdx,
    lcms_uvvis_peak_jcamp_from_df,
)


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

    def tf_jcamp(self):
        if not self.data:
            return None
        
        if self._should_refresh_jcamp() and not self._peaks_applied:
            jdx_file = self._find_peak_or_edit_file()
            if not jdx_file and self.data:
                jdx_file = self.data[0]
            
            if jdx_file:
                lc_df = lcms_df_from_peak_jdx(jdx_file.name)
                if lc_df is None or lc_df.empty:
                    lc_df = lcms_df_from_uvvis_jdx(jdx_file.name)
                if (lc_df is None or lc_df.empty) and self.data:
                    for candidate in self.data:
                        if candidate is jdx_file:
                            continue
                        maybe_df = lcms_df_from_uvvis_jdx(candidate.name)
                        if maybe_df is not None and not maybe_df.empty:
                            lc_df = maybe_df
                            break
                if lc_df is not None and not lc_df.empty:
                    title = self.params.get('fname', 'lc_ms_spectrum')
                    if title:
                        title = title.replace('.jdx', '').replace('.edit', '').replace('.peak', '')
                    else:
                        title = 'lc_ms_spectrum'
                    new_file = lcms_uvvis_peak_jcamp_from_df(lc_df, title, self.params)
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
