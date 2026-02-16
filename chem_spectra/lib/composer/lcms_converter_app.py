import tempfile
from typing import List, Optional, Dict

from chem_spectra.lib.external.chemotion_converter_lcms import (
    lcms_preview_image_from_jdx_files,
    lcms_df_from_peak_jdx,
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

    def tf_img(self):
        if self.params and self.params.get('peaks_str') and not self._peaks_applied:
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
        
        if self.params and self.params.get('peaks_str') and not self._peaks_applied:
            jdx_file = self._find_peak_or_edit_file()
            if not jdx_file and self.data:
                jdx_file = self.data[0]
            
            if jdx_file:
                lc_df = lcms_df_from_peak_jdx(jdx_file.name)
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
