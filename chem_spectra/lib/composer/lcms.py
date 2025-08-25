import tempfile  # noqa: E402

from chem_spectra.lib.composer.base import BaseComposer  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.path as mpath  # noqa: E402
import numpy as np  # noqa: E402


TEXT_SPECTRUM_ORIG = '$$ === CHEMSPECTRA SPECTRUM ORIG ===\n'
TEXT_HEADER_INTEGRATION = '$$ === CHEMSPECTRA INTEGRALS ===\n'
TEXT_PEAK_TABLE_EDIT = '$$ === CHEMSPECTRA PEAK TABLE EDIT ===\n'

class LCMSComposer(BaseComposer):
    def __init__(self, core, lcms_peaks=None, lcms_integrals=None):
        super().__init__(core)
        self.title = core.fname
        if lcms_peaks:
            self.core.edit_peaks = lcms_peaks
        if lcms_integrals:
            self.core.edit_integrals = lcms_integrals
        self.files = self.__compose()
        self.data = self.files

    def __compose(self):
        tic_positive_data, tic_negative_data, uvvis_data, spectra_positive_data, spectra_negative_data = self.core.data
        tic_positive = self.__gen_tic(tic_data=tic_positive_data)
        tic_negative = self.__gen_tic(tic_data=tic_negative_data, is_negative=True)
        tic_positive_jcamp, tic_negative_jcamp = self.tf_jcamp(tic_positive), self.tf_jcamp(tic_negative)

        uvvis_basic = self.__gen_uvvis(data=uvvis_data, include_edits=False)
        uvvis_enriched = self.__gen_uvvis(data=uvvis_data, include_edits=True)
        uvvis_jcamp = self.tf_jcamp(uvvis_basic)
        uvvis_jcamp_enriched = self.tf_jcamp(uvvis_enriched)

        mz_positive = self.__gen_mz_spectra(data=spectra_positive_data)
        mz_positive_jcamp = self.tf_jcamp(mz_positive)

        mz_negative = self.__gen_mz_spectra(data=spectra_negative_data, is_negative=True)
        mz_negative_jcamp = self.tf_jcamp(mz_negative)
        return [uvvis_jcamp_enriched, uvvis_jcamp, tic_positive_jcamp, tic_negative_jcamp, mz_positive_jcamp, mz_negative_jcamp]

    def __gen_tic(self, tic_data, is_negative=False):
        time  = tic_data.get('time', [])
        inten = tic_data.get('Intensity', [])

        if not time:
            category = ('##$CSCATEGORY=TIC NEGATIVE SPECTRUM\n'
                        if is_negative else
                        '##$CSCATEGORY=TIC POSITIVE SPECTRUM\n')
            return [
                '\n', TEXT_SPECTRUM_ORIG,
                f'##TITLE={self.title}\n',
                '##JCAMP-DX=5.00\n',
                '##DATA TYPE=LC/MS\n',
                '##DATA CLASS=XYDATA\n',
                '##ORIGIN=\n', '##OWNER=\n', '##SPECTROMETER/DATA SYSTEM=\n',
                category,
                '##XUNITS=time\n', '##YUNITS=Intensity\n',
                '##XFACTOR=1\n',  '##YFACTOR=1\n',
                '##NPOINTS=0\n',
                '##XYDATA= (XY..XY)\n',
                *self.__gen_ending()
            ]
        max_t, min_t = max(time), min(time)
        max_i, min_i = max(inten), min(inten)
        category = ('##$CSCATEGORY=TIC NEGATIVE SPECTRUM\n'
                    if is_negative else
                    '##$CSCATEGORY=TIC POSITIVE SPECTRUM\n')
        content = [
            '\n', TEXT_SPECTRUM_ORIG,
            f'##TITLE={self.title}\n',
            '##JCAMP-DX=5.00\n',
            '##DATA TYPE=LC/MS\n',
            '##DATA CLASS=XYDATA\n',
            '##ORIGIN=\n', '##OWNER=\n',
            '##SPECTROMETER/DATA SYSTEM=\n',
            category,
            '##XUNITS=time\n', '##YUNITS=Intensity\n',
            '##XFACTOR=1\n',  '##YFACTOR=1\n',
            f'##FIRSTX={time[0]}\n',   f'##LASTX={time[-1]}\n',
            f'##MAXX={max_t}\n',       f'##MAXY={max_i}\n',
            f'##MINX={min_t}\n',       f'##MINY={min_i}\n',
            f'##NPOINTS={len(time)}\n',
            '##XYDATA= (XY..XY)\n',
        ]
        for t, i in zip(time, inten):
            content.append(f'{t}, {i}\n')
        content.extend(self.__gen_ending())
        return content

    def __gen_uvvis(self, data, include_edits=False):
        def _norm_wl(wl):
            try:
                return str(int(round(float(wl))))
            except Exception:
                return str(wl).strip()

        if include_edits:
            raw_user_peaks = getattr(self.core, "edit_peaks", {}) or {}
            user_peaks_by_wl = {_norm_wl(k): v for k, v in raw_user_peaks.items()}

            raw_user_integrals = getattr(self.core, "edit_integrals", {}) or {}
            user_integrals_by_wl = {_norm_wl(k): v for k, v in raw_user_integrals.items()}
        else:
            user_peaks_by_wl = {}
            user_integrals_by_wl = {}

        final_content = []

        for wl_key, value in data.items():
            wl_norm = _norm_wl(wl_key)
            xs, ys = value["RetentionTime"], value["DetectorSignal"]
            
            peaks_to_insert = user_peaks_by_wl.get(wl_norm, []) if include_edits else []
            integrals_to_insert = user_integrals_by_wl.get(wl_norm, []) if include_edits else []

            page_block = [
                '\n',
                TEXT_SPECTRUM_ORIG,
                f"##TITLE={self.title}\n",
                "##JCAMP-DX=5.00\n",
                "##DATA TYPE=LC/MS\n",
                "##DATA CLASS= NTUPLES\n",
                "##ORIGIN=\n",
                "##OWNER=\n",
                "##SPECTROMETER/DATA SYSTEM=\n",
                "##$CSCATEGORY=UVVIS SPECTRUM\n",
                "##VAR_NAME= RETENTION TIME, DETECTOR SIGNAL, WAVELENGTH\n",
                "##SYMBOL= X, Y, T\n",
                "##VAR_TYPE= INDEPENDENT, DEPENDENT, INDEPENDENT\n",
                "##VAR_FORM= AFFN, AFFN, AFFN\n",
                "##VAR_DIM= , , 3\n",
                "##UNITS= RETENTION TIME, DETECTOR SIGNAL, WAVELENGTH\n",
                "##FIRST= , , 1\n",
                f"##PAGE={wl_key}\n",
                f"##NPOINTS={len(xs)}\n",
                "##DATA TABLE= (XY..XY), PEAKS\n",
            ]
            for x, y in zip(xs, ys):
                page_block.append(f"{x}, {y};\n")
            final_content.extend(page_block)

            if include_edits:
                max_x, min_x = (max(xs), min(xs)) if xs else (0, 0)
                max_y, min_y = (max(ys), min(ys)) if ys else (0, 0)
                peak_tbl = [
                    TEXT_PEAK_TABLE_EDIT,
                    f"##TITLE={self.title}\n",
                    "##JCAMP-DX=5.00\n",
                    "##DATA TYPE=UVVISPEAKTABLE\n",
                    "##DATA CLASS=PEAKTABLE\n",
                    "##$CSCATEGORY=EDIT_PEAK\n",
                    "##$CSTHRESHOLD=0.05\n",
                    f"##MAXX={max_x}\n",
                    f"##MAXY={max_y}\n",
                    f"##MINX={min_x}\n",
                    f"##MINY={min_y}\n",
                    "##$CSSOLVENTNAME=\n",
                    "##$CSSOLVENTVALUE=0\n",
                    "##$CSSOLVENTX=0\n",
                    f"##NPOINTS={len(peaks_to_insert)}\n",
                    "##PEAKTABLE= (XY..XY)\n",
                ]
                for p in peaks_to_insert:
                    peak_tbl.append(f"{p['x']}, {p['y']};\n")
                if integrals_to_insert:
                    integrals_str = ' '.join(
                        [f'({",".join(map(str, i))})' for i in integrals_to_insert]
                    )
                    peak_tbl.append(TEXT_HEADER_INTEGRATION)
                    peak_tbl.append(f"##$OBSERVEDINTEGRALS= {integrals_str}\n")

                final_content.extend(peak_tbl)
            
            final_content.append("##END=\n")

        return final_content

    def __gen_mz_spectra(self, data, is_negative=False):
        category = '##$CSCATEGORY=MZ NEGATIVE SPECTRUM\n' if is_negative else '##$CSCATEGORY=MZ POSITIVE SPECTRUM\n'
        content = [
            '\n',
            TEXT_SPECTRUM_ORIG,
            f'##TITLE={self.title}\n',
            '##JCAMP-DX=5.00\n',
            '##DATA TYPE=LC/MS\n',
            '##DATA CLASS= NTUPLES\n',
            '##ORIGIN=\n',
            '##OWNER=\n',
            '##SPECTROMETER/DATA SYSTEM=\n',
            category,
            '##VAR_NAME= M/Z, INTENSITY, RETENTION TIME\n',
            '##SYMBOL= X, Y, T\n',
            '##VAR_TYPE= INDEPENDENT, DEPENDENT, INDEPENDENT\n',
            '##VAR_FORM= AFFN, AFFN, AFFN\n',
            '##VAR_DIM= , , 3\n',
            '##UNITS= M/Z, INTENSITY, RETENTION TIME\n',
            '##FIRST= , , 1\n',
        ]
        
        msspcs_blocks = []
        for time, value in data.items():
            xs, ys = value['mz'], value['intensities']
            msspc = [
                f'##PAGE={time}\n',
                f'##NPOINTS={len(xs)}\n',
                '##DATA TABLE= (XY..XY), PEAKS\n',
            ]
            for x, y in zip(xs, ys):
                msspc.append(f'{x}, {y};\n')
            msspcs_blocks.append(''.join(msspc))

        msspcs = '\n'.join(msspcs_blocks)
        content.append(msspcs)
        content.extend(self.__gen_ending())

        return content
        
    def __gen_ending(self):
        return [
            '\n##END=\n',
            '\n'
        ]
    
    def tf_jcamp(self, data):
        meta = ''.join(data)
        tf = tempfile.NamedTemporaryFile(suffix='.jdx')
        tf.write(bytes(meta, 'UTF-8'))
        tf.seek(0)
        return tf

    def tf_img(self):
        plt.rcParams['figure.figsize'] = [16, 9]
        plt.rcParams['figure.dpi'] = 200
        plt.rcParams['font.size'] = 14

        try:
            _, _, uvvis_data, _, _ = self.core.data
        except Exception:
            uvvis_data = {}

        if uvvis_data:
            keys = list(uvvis_data.keys())

            def _to_float_or_none(k):
                try:
                    return float(k)
                except Exception:
                    return None

            numeric_pairs = [(k, _to_float_or_none(k)) for k in keys]
            numeric_only = [p for p in numeric_pairs if p[1] is not None]
            wl_key = min(numeric_only, key=lambda p: p[1])[0] if numeric_only else keys[0]

            wl_entry = uvvis_data.get(wl_key, {"RetentionTime": [], "DetectorSignal": []})
            xs = wl_entry.get('RetentionTime') or []
            ys = wl_entry.get('DetectorSignal') or []

            xs = np.asarray(xs, float) if len(xs) else np.asarray([])
            ys = np.asarray(ys, float) if len(ys) else np.asarray([])

            plt.plot(xs, ys)
            plt.xlabel('X (Retention Time)', fontsize=18)
            plt.ylabel('Y (Detector Signal)', fontsize=18)
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
            y_boundary_min = y_min - h * 0.05
            y_boundary_max = y_max + h * 0.15

            def _norm_wl(wl):
                try:
                    return str(int(round(float(wl))))
                except Exception:
                    return str(wl).strip()

            raw_user_peaks = getattr(self.core, "edit_peaks", {}) or {}
            user_peaks_by_wl = {_norm_wl(k): v for k, v in raw_user_peaks.items()}
            wl_norm = _norm_wl(wl_key)
            peaks_to_plot = user_peaks_by_wl.get(wl_norm) or user_peaks_by_wl.get('default') or []

            def _detect_time_scale(peaks_list, series_x):
                try:
                    xs_max = float(np.max(series_x)) if series_x.size else 0.0
                    px_vals = [float(p.get('x')) for p in peaks_list if p and p.get('x') is not None]
                    px_max = max(px_vals) if px_vals else 0.0
                    if xs_max > 300 and 0 < px_max <= 60:
                        return 60.0
                except Exception:
                    pass
                return 1.0

            scale_factor = _detect_time_scale(peaks_to_plot, xs)

            path_data = [(mpath.Path.MOVETO, (0, 0)), (mpath.Path.LINETO, (0, 1))]
            codes, verts = zip(*path_data)
            marker = mpath.Path(verts, codes)

            x_peaks, y_peaks = [], []
            if isinstance(peaks_to_plot, list) and peaks_to_plot and xs.size and ys.size:
                def _nearest_y(x_target: float) -> float:
                    idx = int(np.argmin(np.abs(xs - x_target)))
                    return float(ys[idx])

                for p in peaks_to_plot:
                    try:
                        x_val = float(p.get('x')) * scale_factor
                        y_at_x = _nearest_y(x_val)
                        x_peaks.append(x_val)
                        y_peaks.append(y_at_x)
                    except Exception:
                        continue

            if x_peaks:
                y_offset = 0.02 * h
                y_peaks_plot = [y + y_offset for y in y_peaks]
                plt.plot(x_peaks, y_peaks_plot, 'r', ls='', marker=marker, markersize=50, zorder=3)
                y_boundary_max = max(y_boundary_max, max(y_peaks_plot) + 0.05 * h)

            raw_user_integrals = getattr(self.core, "edit_integrals", {}) or {}
            user_integrals_by_wl = {_norm_wl(k): v for k, v in raw_user_integrals.items()}
            integrals_to_plot = user_integrals_by_wl.get(wl_norm) or []

            if isinstance(integrals_to_plot, list) and integrals_to_plot and xs.size and ys.size:
                def _find_idx(x_target: float) -> int:
                    return int(np.argmin(np.abs(xs - x_target)))

                for itg in integrals_to_plot:
                    if not (isinstance(itg, (list, tuple)) and len(itg) >= 2):
                        continue
                    try:
                        xL = float(itg[0]) * scale_factor
                        xU = float(itg[1]) * scale_factor
                    except Exception:
                        continue

                    iL, iU = _find_idx(xL), _find_idx(xU)
                    if iL == iU:
                        continue
                    if iL > iU:
                        iL, iU = iU, iL

                    cxs = xs[iL:iU + 1]
                    cys = ys[iL:iU + 1]
                    if cxs.size < 2:
                        continue

                    yL, yU = float(cys[0]), float(cys[-1])
                    slope = (yU - yL) / (float(cxs[-1]) - float(cxs[0]) + 1e-12)
                    baseline = slope * (cxs - cxs[0]) + yL

                    plt.fill_between(cxs, y1=cys, y2=baseline, alpha=0.20, color='#FF0000', zorder=1)

            plt.ylim(y_boundary_min, y_boundary_max)

        tf_img = tempfile.NamedTemporaryFile(suffix='.png')
        plt.savefig(tf_img, format='png')
        tf_img.seek(0)
        plt.clf()
        plt.cla()
        plt.close()
        return tf_img
