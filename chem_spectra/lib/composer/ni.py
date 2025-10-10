import json
from operator import truediv
import uuid
import matplotlib
matplotlib.use('Agg')

import tempfile  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.path as mpath  # noqa: E402
import re  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import ticker  # noqa: E402
from matplotlib.patches import FancyArrowPatch
import csv

from chem_spectra.lib.composer.base import (  # noqa: E402, F401
    extrac_dic, calc_npoints, BaseComposer
)
from chem_spectra.lib.shared.calc import (  # noqa: E402
    calc_mpy_center, calc_ks, get_curve_endpoint, cal_slope, cal_xyIntegration,
)


TEXT_SPECTRUM_ORIG = '$$ === CHEMSPECTRA SPECTRUM ORIG ===\n'
TEXT_INTEGRATION = '$$ === CHEMSPECTRA INTEGRATION ===\n'
TEXT_MULTIPLICITY = '$$ === CHEMSPECTRA MULTIPLICITY ===\n'


class NIComposer(BaseComposer):
    def __init__(self, core):
        super().__init__(core)
        self.title = core.fname
        self.__override_cv_density_label()
        self.meta = self.__compose()

    def __override_cv_density_label(self):
        try:
            if not getattr(self.core, 'is_cyclic_volta', False):
                return
            cv_params = self.core.params.get('cyclicvolta', {}) if hasattr(self.core, 'params') else {}
            if isinstance(cv_params, str):
                cv_params = json.loads(cv_params)
            use_density = bool(cv_params.get('useCurrentDensity', False))
            if not use_density:
                return

            area_unit = cv_params.get('areaUnit') or 'cm²'
            axes_units = self.core.params.get('axesUnits') if hasattr(self.core, 'params') else None
            base_unit = 'A'
            if isinstance(axes_units, dict):
                selected_y = str(axes_units.get('yUnit') or '')
                base_unit = 'mA' if ('mA' in selected_y or 'mA' in selected_y) else 'A'

            self.core.label['y'] = f"Current density in {base_unit}/{area_unit}"
        except Exception:
            pass

    def __header_base(self):
        return [
            '\n',
            TEXT_SPECTRUM_ORIG,
            '##TITLE={}\n'.format(self.title),
            '##JCAMP-DX=5.00\n',
            '##DATA TYPE={}\n'.format(self.core.datatype),
            '##DATA CLASS=XYDATA\n',
            '##$CSCATEGORY=SPECTRUM\n',
            '##ORIGIN={}\n'.format(extrac_dic(self.core, 'ORIGIN')),
            '##OWNER={}\n'.format(extrac_dic(self.core, 'OWNER')),
        ]

    def __calc_nucleus_by_boundary(self):
        delta = self.core.boundary['x']['max'] - self.core.boundary['x']['min']
        nucleus = '^1H' if abs(delta) < 40.0 else '^13C'
        return nucleus

    def __get_nucleus(self):
        nuc_orig = extrac_dic(self.core, '.OBSERVENUCLEUS')
        nuc_modf = re.sub('[^A-Za-z0-9]+', '', nuc_orig).lower()
        is_valid = ('13c' in nuc_modf) or ('1h' in nuc_modf) or ('19f' in nuc_modf) or ('31p' in nuc_modf) or ('15n' in nuc_modf) or ('29si' in nuc_modf)   # noqa: E501
        nucleus = nuc_orig if is_valid else self.__calc_nucleus_by_boundary()
        return nucleus

    def __header_nmr(self):
        header_lines = []

        # Append each line only if the extracted value is not empty
        observe_frequency = extrac_dic(self.core, '.OBSERVEFREQUENCY')
        if observe_frequency:
            header_lines.append('##.OBSERVE FREQUENCY={}\n'.format(observe_frequency))

        observe_nucleus = self.__get_nucleus()
        if observe_nucleus:
            header_lines.append('##.OBSERVE NUCLEUS={}\n'.format(observe_nucleus))

        spectrometer_data_system = extrac_dic(self.core, 'SPECTROMETER/DATASYSTEM')
        if spectrometer_data_system:
            header_lines.append('##SPECTROMETER/DATA SYSTEM={}\n'.format(spectrometer_data_system))

        shift_reference = extrac_dic(self.core, '.SHIFTREFERENCE')
        if shift_reference:
            header_lines.append('##.SHIFT REFERENCE={}\n'.format(shift_reference))

        solvent_name = extrac_dic(self.core, '.SOLVENTNAME')
        if solvent_name:
            header_lines.append('##.SOLVENT NAME={}\n'.format(solvent_name))

        pulse_sequence = extrac_dic(self.core, '.PULSESEQUENCE')
        if pulse_sequence:
            header_lines.append('##.PULSE SEQUENCE={}\n'.format(pulse_sequence))

        return header_lines

    def __header_params(self):
        return [
            '##XUNITS={}\n'.format(self.core.label['x']),
            '##YUNITS={}\n'.format(self.core.label['y']),
            '##XFACTOR={}\n'.format(self.core.factor['x']),
            '##YFACTOR={}\n'.format(self.core.factor['y']),
            '##FIRSTX={}\n'.format(self.core.first_x),
            '##LASTX={}\n'.format(self.core.last_x),
            '##MAXX={}\n'.format(self.core.boundary['x']['max']),
            '##MAXY={}\n'.format(self.core.boundary['y']['max']),
            '##MINX={}\n'.format(self.core.boundary['x']['min']),
            '##MINY={}\n'.format(self.core.boundary['y']['min'])
        ]

    def __gen_headers_spectrum_orig(self):
        if self.core.is_em_wave or self.core.non_nmr:
            return self.__header_base() + self.__header_params()
        else:
            return self.__header_base() + \
                self.__header_nmr() + self.__header_params()

    def __gen_headers_im(self):
        return [
            '$$ === CHEMSPECTRA INTEGRALS AND MULTIPLETS ===\n',
        ]

    def __gen_headers_integration(self):
        return [
            '##$OBSERVEDINTEGRALS= (X Y Z)\n',
        ]

    def __gen_headers_mpy_integ(self):
        return [
            '##$OBSERVEDMULTIPLETS=\n',
        ]

    def __gen_headers_mpy_peaks(self):
        return [
            '##$OBSERVEDMULTIPLETSPEAKS=\n',
        ]

    def __gen_header_simulation(self):
        return [
            '$$ === CHEMSPECTRA SIMULATION ===\n',
            '##$CSSIMULATIONPEAKS=\n',
        ]

    def __gen_header_cyclic_voltammetry(self):
        return [
            '$$ === CHEMSPECTRA CYCLIC VOLTAMMETRY ===\n',
        ]

    def __gen_header_sec(self):
        core_dic = self.core.dic
        sec_data_key = ['MN', 'MW', 'MP', 'D']
        result = []
        for key in sec_data_key:
            dic_value = core_dic.get(key, [])
            key_str = f"##{key}={dic_value[0]}\n" if len(dic_value) > 0 else f'##{key}=\n'
            result.append(key_str)
        return result

    def __gen_header_user_input_meta_data(self):
        if self.core.is_dsc:
            dsc_meta_data = self.core.params.get('dsc_meta_data', None)
            melting_point, tg_value = '', ''
            if dsc_meta_data is not None:
                melting_point = dsc_meta_data.get('meltingPoint', '')
                tg_value = dsc_meta_data.get('tg', '')
            else:
                melting_point_arr = self.core.dic.get('MELTINGPOINT', [''])
                tg_value_arr = self.core.dic.get('TG', [''])
                melting_point = melting_point_arr[0]
                tg_value = tg_value_arr[0]
            return [
                f'##MELTINGPOINT={melting_point}\n',
                f'##TG={tg_value}\n'
            ]
        return []

    def __get_xy_of_peak(self, peak):
        if peak is None:
            return '', ''
        x = peak['x'] or ''
        y = peak['y'] or ''
        return x, y

    def __gen_cyclic_voltammetry_medadata(self):
        meta = []
        scan_rate = self.core.dic.get('SCANRATE', [0.1])[0]
        x_values = self.core.xs
        spectrum_direction = ''
        if len(x_values) > 2:
            spectrum_direction = 'NEGATIVE' if x_values[0] > x_values[1] else 'POSITIVE'

        cv_params = {}
        if hasattr(self.core, 'params'):
            cv_params = self.core.params.get('cyclicvolta', {}) or {}
            if isinstance(cv_params, str):
                cv_params = json.loads(cv_params)

        area_value = cv_params.get('areaValue', '')
        area_unit = cv_params.get('areaUnit', '')
        use_current_density = cv_params.get('useCurrentDensity', False)

        axis_length_unit = ''
        if isinstance(area_unit, str):
            if 'cm' in area_unit:
                axis_length_unit = 'cm'
            elif 'mm' in area_unit:
                axis_length_unit = 'mm'
        
        meta.append(f"##$CSSCANRATE={scan_rate}\n")
        meta.append(f"##$CSSPECTRUMDIRECTION={spectrum_direction}\n")
        meta.append(f"##$CSWEAREAVALUE={area_value}\n")
        meta.append(f"##$CSWEAREAUNIT={area_unit}\n")
        meta.append(f"##$CSCURRENTMODE={'DENSITY' if use_current_density else 'CURRENT'}\n")
        return meta

    def __gen_cyclic_voltammetry_data_peaks(self):
        content = ['##$CSCYCLICVOLTAMMETRYDATA=\n']
        if self.core.is_cyclic_volta:
            listMaxMinPeaks = self.core.max_min_peaks_table
            cyclicvolta_data = self.core.params['cyclicvolta']
            current_jcamp_idx = self.core.params['jcamp_idx']
            if self.core.params['list_max_min_peaks'] is not None:
                listMaxMinPeaks = self.core.params['list_max_min_peaks']

            x_max = np.max(self.core.xs)
            arr_max_x_indices = np.where(self.core.xs == x_max)
            y_pecker = 0
            if (arr_max_x_indices[0][0]):
                idx = arr_max_x_indices[0][0]
                y_pecker = self.core.ys[idx]

            for peak in listMaxMinPeaks:
                max_peak, min_peak = None, None
                if 'max' in peak:
                    max_peak = peak['max']
                if 'min' in peak:
                    min_peak = peak['min']
                x_max_peak, y_max_peak = self.__get_xy_of_peak(max_peak)
                x_min_peak, y_min_peak = self.__get_xy_of_peak(min_peak)

                x_pecker = ''
                if 'pecker' in peak and peak['pecker'] is not None:
                    pecker = peak['pecker']
                    x_pecker = pecker['x']
                    y_pecker = pecker['y']
                    x_pecker = f"{float(x_pecker)}"

                if (x_max_peak == '' and x_min_peak == '' and x_pecker == ''):
                    # ignore if missing both max and min peak
                    continue

                if (x_max_peak == '' or x_min_peak == ''):
                    delta = ''
                else:
                    delta = abs(x_max_peak - x_min_peak)

                # calculate ratio
                if (y_min_peak == '' or y_max_peak == ''):
                    ratio = ''
                else:
                    first_expr = abs(y_min_peak) / abs(y_max_peak)
                    second_expr = 0.485 * abs(y_pecker) / abs(y_max_peak)
                    ratio = first_expr + second_expr + 0.086
                    if (y_pecker) == 0:
                        y_pecker = ''

                is_ref = peak.get('isRef', False)
                is_ref_in_number = 1 if is_ref else 0
                content.append(
                    '({x_max}, {y_max}, {x_min}, {y_min}, {ratio}, {delta}, {x_pecker}, {y_pecker}, {is_ref})\n'.format(x_max=x_max_peak, y_max=y_max_peak, x_min=x_min_peak, y_min=y_min_peak, ratio=ratio, delta=delta, x_pecker=x_pecker, y_pecker=y_pecker, is_ref=is_ref_in_number)  # noqa: E501
                )

        return content

    def __compose(self):
        meta = []
        meta.extend(self.gen_headers_root())

        meta.extend(self.__gen_headers_spectrum_orig())
        if self.core.is_sec:
            meta.extend(self.__gen_header_sec())
        meta.extend(self.__gen_header_user_input_meta_data())
        meta.extend(self.gen_spectrum_orig())
        meta.extend(self.__gen_headers_im())
        meta.extend(self.__gen_headers_integration())
        meta.extend(self.gen_integration_info())
        meta.extend(self.__gen_headers_mpy_integ())
        meta.extend(self.gen_mpy_integ_info())
        meta.extend(self.__gen_headers_mpy_peaks())
        meta.extend(self.gen_mpy_peaks_info())
        meta.extend(self.__gen_header_simulation())
        meta.extend(self.gen_simulation_info())
        if self.core.is_cyclic_volta:
            meta.extend(self.__gen_header_cyclic_voltammetry())
            meta.extend(self.__gen_cyclic_voltammetry_medadata())
            meta.extend(self.__gen_cyclic_voltammetry_data_peaks())
        meta.extend(self.gen_ending())

        meta.extend(self.gen_headers_peaktable_edit())
        meta.extend(self.gen_edit_peaktable())
        meta.extend(self.gen_ending())

        meta.extend(self.gen_headers_peaktable_auto())
        meta.extend(self.gen_auto_peaktable())
        meta.extend(self.gen_ending())

        meta.extend(self.generate_original_metadata())

        meta.extend(self.gen_ending())
        return meta

    def __plt_nbins(self):
        return 20

    def __fakto(self):
        typ = self.core.typ
        if 'INFRARED' == typ:
            return -1
        return 1

    def tf_img(self):
        plt.rcParams['figure.figsize'] = [16, 9]
        plt.rcParams['figure.dpi'] = 200
        plt.rcParams['font.size'] = 14

        # PLOT data
        plt.plot(self.core.xs, self.core.ys)
        x_max, x_min = self.core.boundary['x']['max'], self.core.boundary['x']['min']   # noqa: E501

        xlim_left, xlim_right = [x_min, x_max] if (self.core.is_tga or self.core.is_gc or self.core.is_uv_vis or self.core.is_hplc_uv_vis or self.core.is_xrd or self.core.is_cyclic_volta or self.core.is_sec or self.core.is_cds or self.core.is_aif or self.core.is_emissions or self.core.is_dls_acf or self.core.is_dls_intensity) else [x_max, x_min]    # noqa: E501
        plt.xlim(xlim_left, xlim_right)
        y_max, y_min = np.max(self.core.ys), np.min(self.core.ys)
        h = y_max - y_min
        w = x_max - x_min
        y_boundary_min = y_min - h * 0.2
        y_boundary_max = y_max + h * 0.5

        # PLOT peaks
        faktor = self.__fakto()
        path_data = [
            (mpath.Path.MOVETO, (0, faktor * 5)),
            (mpath.Path.LINETO, (0, faktor * 20)),
        ]
        codes, verts = zip(*path_data)
        marker = mpath.Path(verts, codes)

        circle = mpath.Path.unit_circle()
        cirle_verts = np.concatenate([circle.vertices, verts])
        cirle_codes = np.concatenate([circle.codes, codes])
        cut_star_marker = mpath.Path(cirle_verts, cirle_codes)

        x_peaks = []
        y_peaks = []
        if self.core.edit_peaks:
            x_peaks = self.core.edit_peaks['x']
            y_peaks = self.core.edit_peaks['y']
        elif self.core.auto_peaks:
            x_peaks = self.core.auto_peaks['x']
            y_peaks = self.core.auto_peaks['y']

        x_peckers = []
        y_peckers = []
        x_peaks_ref, y_peaks_ref = [], []
        if self.core.is_cyclic_volta:
            x_peaks = []
            y_peaks = []
            listMaxMinPeaks = []
            if self.core.params['list_max_min_peaks'] is not None:
                listMaxMinPeaks = self.core.params['list_max_min_peaks']

            for peak in listMaxMinPeaks:
                max_peak, min_peak = None, None
                if 'max' in peak:
                    max_peak = peak['max']
                if 'min' in peak:
                    min_peak = peak['min']
                x_max_peak, y_max_peak = self.__get_xy_of_peak(max_peak)
                x_min_peak, y_min_peak = self.__get_xy_of_peak(min_peak)

                if (x_max_peak == '' and x_min_peak == ''):
                    # ignore if missing both max and min peak
                    continue

                if (x_max_peak == '' and y_max_peak == ''):
                    x_peaks.extend([x_min_peak])
                    y_peaks.extend([y_min_peak])
                elif (x_min_peak == '' and y_min_peak == ''):
                    x_peaks.extend([x_max_peak])
                    y_peaks.extend([y_max_peak])
                else:
                    is_ref = peak.get('isRef', False) if 'isRef' in peak else False
                    if is_ref:
                        x_peaks_ref.extend([x_max_peak, x_min_peak])
                        y_peaks_ref.extend([y_max_peak, y_min_peak])
                    else:
                        x_peaks.extend([x_max_peak, x_min_peak])
                        y_peaks.extend([y_max_peak, y_min_peak])

                if 'pecker' in peak and peak['pecker'] is not None:
                    pecker = peak['pecker']
                    x_pecker, y_pecker = pecker['x'], pecker['y']
                    x_peckers.append(x_pecker)
                    y_peckers.append(y_pecker)

            # display x value of peak for cyclic voltammetry
            for i in range(len(x_peaks)):
                x_pos = x_peaks[i]
                y_pos = y_peaks[i] + h * 0.1
                x_float = '{:.2e}'.format(x_pos)
                y_float = '{:.2e}'.format(y_peaks[i])
                peak_label = 'x: {x}\ny: {y}'.format(x=x_float, y=y_float)
                plt.text(x_pos, y_pos, peak_label)

            # display x value of ref peak for cyclic voltammetry
            for i in range(len(x_peaks_ref)):
                x_pos = x_peaks_ref[i]
                y_pos = y_peaks_ref[i] + h * 0.1
                x_float = '{:.2e}'.format(x_pos)
                y_float = '{:.2e}'.format(y_peaks_ref[i])
                peak_label = 'x: {x}\ny: {y}'.format(x=x_float, y=y_float)
                plt.text(x_pos, y_pos, peak_label)

        plt.plot(
            x_peaks,
            y_peaks,
            'r',
            ls='',
            marker=marker,
            markersize=50,
        )

        plt.plot(
            x_peckers,
            y_peckers,
            'g',
            ls='',
            marker=marker,
            markersize=50,
        )

        plt.plot(
            x_peaks_ref,
            y_peaks_ref,
            'r',
            ls='',
            marker=cut_star_marker,
            markersize=50,
        )

        # ----- Calculate integration -----
        refShift, refArea = self.refShift, self.refArea
        if (len(self.all_itgs) == 0 and len(self.core.itg_table) > 0 and not self.core.params['integration'].get('edited') and ('originStack' not in self.core.params['integration'])):
            core_itg_table = self.core.itg_table[0]
            itg_table = core_itg_table.split('\n')
            for itg in itg_table:
                clear_itg = itg.replace('(', '')
                clear_itg = clear_itg.replace(')', '')
                split_itg = clear_itg.split(',')
                if (len(split_itg) > 2):
                    xLStr = split_itg[0].strip()
                    xUStr = split_itg[1].strip()
                    areaStr = split_itg[2].strip()
                    self.all_itgs.append({'xL': float(xLStr), 'xU': float(xUStr), 'area': float(areaStr)})    # noqa: E501
        

        # ----- Calculate multiplicity -----
        if (len(self.mpys) == 0 and len(self.core.mpy_itg_table) > 0 and not self.core.params['integration'].get('edited') and ('originStack' not in self.core.params['integration'])):
            core_mpy_pks_table = self.core.mpy_pks_table[0]
            mpy_pks_table = core_mpy_pks_table.split('\n')
            tmp_dic_mpy_peaks = {}
            for peak in mpy_pks_table:
                clear_peak = peak.replace('(', '')
                clear_peak = clear_peak.replace(')', '')
                split_peak = clear_peak.split(',')
                idx_peakStr = split_peak[0].strip()
                xStr = split_peak[1].strip()
                yStr = split_peak[2].strip()
                if idx_peakStr not in tmp_dic_mpy_peaks:
                    tmp_dic_mpy_peaks[idx_peakStr] = []
                
                tmp_dic_mpy_peaks[idx_peakStr].append({'x': float(xStr), 'y': float(yStr)})

            core_mpy_itg_table = self.core.mpy_itg_table[0]
            mpy_itg_table = core_mpy_itg_table.split('\n')
            for mpy in mpy_itg_table:
                clear_mpy = mpy.replace('(', '')
                clear_mpy = clear_mpy.replace(')', '')
                split_mpy = clear_mpy.split(',')
                mpy_item = { 'mpyType': '', 'xExtent': {'xL': 0.0, 'xU': 0.0}, 'yExtent': {'yL': 0.0, 'yU': 0.0}, 'peaks': [], 'area': 1.0 }
                if (len(split_mpy) > 7):
                    idxStr = split_mpy[0].strip()
                    xLStr = split_mpy[1].strip()
                    xUStr = split_mpy[2].strip()
                    mpy_item['xExtent']['xL'] = float(xLStr) + refShift
                    mpy_item['xExtent']['xU'] = float(xUStr) + refShift
                    yLStr = split_mpy[3].strip()
                    yUStr = split_mpy[4].strip()
                    mpy_item['yExtent']['yL'] = float(yLStr) + refShift
                    mpy_item['yExtent']['yU'] = float(yUStr) + refShift
                    typeStr = split_mpy[6].strip()
                    mpy_item['mpyType'] = typeStr
                    mpy_item['peaks'] = tmp_dic_mpy_peaks[idxStr]
                    self.mpys.append(mpy_item)

        # ----- PLOT integration -----
        itg_h = y_max + h * 0.6
        itg_h = itg_h + itg_h * 0.1
        itg_value_position_y = y_min - h * 0.25
        if (len(self.mpys) == 0):
            itg_value_position_y =  y_min - h * 0.05
        y_boundary_max = self.__draw_integrals(plt, refShift, refArea, y_max, h, itg_value_position_y, itg_h)
        y_boundary_min = itg_value_position_y - h * 0.1
        
        # ----- PLOT multiplicity -----
        mpy_h = y_min - h * 0.03
        for mpy in self.mpys:
            xL, xU, area, typ, peaks = mpy['xExtent']['xL'] - refShift, mpy['xExtent']['xU'] - refShift, mpy['area'] * refArea, mpy['mpyType'], mpy['peaks']    # noqa: E501
            plt.plot([xL, xU], [mpy_h, mpy_h], color='#DA70D6')
            plt.plot([xL, xL], [mpy_h + h * 0.01, mpy_h - h * 0.01], color='#DA70D6')   # noqa: E501
            plt.plot([xU, xU], [mpy_h + h * 0.01, mpy_h - h * 0.01], color='#DA70D6')   # noqa: E501
            plt.text((xL + xU) / 2, mpy_h - h * 0.01, '{:0.3f} ({})'.format(calc_mpy_center(mpy['peaks'], refShift, mpy['mpyType']), typ), color='#DA70D6', size=7, rotation=90., ha='right', va='top', rotation_mode='anchor')  # noqa: E501
            for p in peaks:
                x = p['x']
                plt.plot([x - refShift, x - refShift], [mpy_h, mpy_h + h * 0.02], color='#DA70D6')  # noqa: E501

        # PLOT label
        if (self.core.is_xrd):
            waveLength = self.core.params['waveLength']
            label = "X ({}), WL={} nm".format(self.core.label['x'], waveLength['value'], waveLength['unit'])    # noqa: E501
            plt.xlabel((label), fontsize=18)
        elif (self.core.is_cyclic_volta):
            plt.xlabel("{}".format(self.core.label['x']), fontsize=18)
        elif (self.core.non_nmr == False):
            plt.xlabel("Chemical shift ({})".format(self.core.label['x'].lower()), fontsize=18)
        else:
            plt.xlabel("X ({})".format(self.core.label['x']), fontsize=18)

        if (self.core.is_cyclic_volta):
            plt.ylabel("{}".format(self.core.label['y']), fontsize=18)
        elif (self.core.non_nmr == False):
            plt.ylabel("Intensity ({})".format(self.core.label['y'].lower()), fontsize=18)
        else:
            plt.ylabel("Y ({})".format(self.core.label['y']), fontsize=18)
        plt.locator_params(nbins=self.__plt_nbins())
        plt.grid(False)

        self.__generate_info_box(plt)

        y_boundary_max = self.__draw_peaks(plt, x_peaks, y_peaks, h, w, y_boundary_max * (1.1 if self.core.is_ir else 1.5))


        plt.ylim(
            y_boundary_min,
            y_boundary_max,
        )

        ax = plt.gca()
        if self.core.is_cyclic_volta:
            cv = (self.core.params.get('cyclicvolta') or {})
            if isinstance(cv, str):
                cv = json.loads(cv)
            if cv.get('useCurrentDensity', False):
                unit = str(cv.get('areaUnit') or 'cm²').lower()
                u2cm2 = {'mm²': 0.01, 'mm2': 0.01, 'cm²': 1.0, 'cm2': 1.0, 'm²': 10000.0, 'm2': 10000.0}
                try:
                    area_val = float(cv.get('areaValue', 1.0))
                except Exception:
                    area_val = 1.0
                divisor = max(area_val * u2cm2.get(unit, 1.0), 1e-30)

                y0, y1 = ax.get_ylim()
                ymax_disp = max(abs(y0), abs(y1)) / divisor
                exp = int(np.floor(np.log10(ymax_disp))) if ymax_disp > 0 else 0
                base = (10.0 ** exp) if exp != 0 else 1.0

                from matplotlib.ticker import FuncFormatter
                ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _:
                    f"{(y/divisor)/base:.3g}"
                ))

                ax.yaxis.get_offset_text().set_visible(False)

                if exp != 0:
                    ax.text(
                        0.0, 1,              
                        r"$\times 10^{%d}$" % exp,
                        transform=ax.transAxes,
                        ha='left', va='bottom',
                        fontsize=14,
                        clip_on=False
                    )
            else:
                fmt = ticker.ScalarFormatter(useMathText=True)
                fmt.set_scientific(True)
                fmt.set_powerlimits((-1, 1))
                ax.yaxis.set_major_formatter(fmt)
        else:
            fmt = ticker.ScalarFormatter(useMathText=True)
            fmt.set_scientific(True)
            fmt.set_powerlimits((-1, 1))
            ax.yaxis.set_major_formatter(fmt)
                    
        # Save
        tf_img = tempfile.NamedTemporaryFile(suffix='.png')
        plt.savefig(tf_img, format='png')
        tf_img.seek(0)
        plt.clf()
        plt.cla()
        return tf_img
    

    def __draw_integrals(self, plt, refShift, refArea, y_max, h, itg_value_position_y, itg_h):
        y_boundary_max = y_max
        for itg in self.all_itgs:
            # integration marker
            xL, xU, area = itg['xL'] - refShift, itg['xU'] - refShift, itg['area'] * refArea    # noqa: E501
            # integration curve
            ks = calc_ks(self.core.ys, y_max, h)
            iL, iU = get_curve_endpoint(self.core.xs, self.core.ys, xL, xU)
            ref = ks[iL]
            cxs = self.core.xs[iL:iU]

            if self.core.is_hplc_uv_vis:
                # fill area under curve
                cys = self.core.ys[iL:iU]
                slope = cal_slope(cxs[0], cys[0], cxs[len(cxs)-1], cys[len(cys)-1])  # noqa: E501
                last_y = cys[0]
                last_x = cxs[0]
                aucys = [last_y]
                for i in range(1, len(cys)):
                    curr_x = cxs[i]
                    curr_y = slope*(curr_x-last_x) + last_y
                    aucys.append(curr_y)
                    last_x = curr_x
                    last_y = curr_y
                plt.fill_between(cxs, y1=cys, y2=aucys, alpha=0.2, color='#FF0000')  # noqa: E501
            else:
                # display integration
                plt.plot([xL, xU], [itg_value_position_y, itg_value_position_y], color='#228B22')
                plt.plot([xL, xL], [itg_value_position_y + h * 0.01, itg_value_position_y - h * 0.01], color='#228B22')   # noqa: E501
                plt.plot([xU, xU], [itg_value_position_y + h * 0.01, itg_value_position_y - h * 0.01], color='#228B22')   # noqa: E501
                plt.text((xL + xU) / 2, itg_value_position_y - h * 0.01, '{:0.2f}'.format(area), color='#228B22', ha='right', rotation_mode='anchor', size=7, rotation=90.)   # noqa: E501
                cys = (ks[iL:iU] - ref) * 1.5 + h * 0.15
                plt.plot(cxs, cys, color='#228B22')
                try:
                    cys_max = np.max(cys)
                    y_boundary_max = max(y_boundary_max, cys_max + h*0.1)
                except:
                    pass

        return y_boundary_max
        

    def __draw_peaks(self, plt, x_peaks, y_peaks, h, w, y_boundary_max):
        if self.core.non_nmr == True or len(x_peaks) == 0:
            return y_boundary_max

        params = self.core.params
        if params['ref_name'] is None or params['peaks_str'] is None or params['peaks_str'] == '':
            return y_boundary_max

        ax = plt.gca()
        differences = np.diff(x_peaks)
        grouping_threshold = np.mean(differences)
        groups_x = []
        groups_y = []
        current_group_x = [x_peaks[0]]
        current_group_y = [y_peaks[0]]
        
        for i in range(1, len(x_peaks)):
            if (x_peaks[i] - current_group_x[-1] < grouping_threshold):
                current_group_x.append(x_peaks[i])
                current_group_y.append(y_peaks[i])
            else:
                groups_x.append(current_group_x)
                current_group_x = [x_peaks[i]]
                groups_y.append(current_group_y)
                current_group_y = [y_peaks[i]]

        groups_x.append(current_group_x)
        groups_y.append(current_group_y)
        max_values = [np.max(items) for items in groups_y]
        x_boundary_min, x_boundary_max = np.min(x_peaks), np.max(x_peaks)

        highest_peak_anotation = 0

        for i in range(len(groups_x)):
            mygroup_x, mygroup_y = groups_x[i], groups_y[i]
            len_group_x = len(mygroup_x)
            if len_group_x > 1:
                max_current_group = max_values[i] + h * 0.25
                if (i > 0):
                    prev_max_group = max_values[i-1] + h * 0.25
                    my_gap = abs(max_current_group - prev_max_group)
                    if my_gap < h*0.1:
                        max_current_group = max_current_group + h * 0.45
                middle_idx = int(len_group_x/2)
                gap_value = np.mean(mygroup_x)
                x_text = 0
                for j in range(len_group_x):
                    x_pos = mygroup_x[j]
                    y_pos = mygroup_y[j] + h * 0.5
                    x_float = '{:.2f}'.format(x_pos)
                    peak_label = '{x}'.format(x=x_float)
                    if j >= middle_idx:
                        x_text = -(w * 0.01) * (middle_idx - j)
                    elif j < middle_idx:
                        x_text = w * 0.01 * (j-middle_idx)

                    ax.add_patch(FancyArrowPatch((x_pos, max_current_group), (x_pos, max_current_group  + h * 0.05), linewidth=0.1))
                    ax.add_patch(FancyArrowPatch((x_pos, max_current_group  + h * 0.05), (gap_value + x_text, max_current_group  + h * 0.11), linewidth=0.1))

                    x_boundary_min = min(gap_value + x_text, x_boundary_min)
                    x_boundary_max = max(gap_value + x_text, x_boundary_max)

                    ax.annotate(peak_label,
                        xy=(gap_value + x_text, max_current_group  + h * 0.11), xycoords='data',
                        xytext=(0, 12), textcoords='offset points',
                        arrowprops=dict(arrowstyle="-", linewidth=0.2),
                        rotation=90, size=6)

                    highest_peak_anotation = max(highest_peak_anotation, max_current_group  + h * 0.25)
            else:
                x_pos = mygroup_x[0]
                y_pos = mygroup_y[0] + h * 0.18
                x_float = '{:.2f}'.format(x_pos)
                peak_label = '{x}'.format(x=x_float)
                ax.annotate(peak_label,
                    xy=(x_pos, y_pos), xycoords='data',
                    xytext=(0, 20), textcoords='offset points',
                    arrowprops=dict(arrowstyle="-", linewidth=0.2),
                    rotation=90, size=6)

                highest_peak_anotation = max(highest_peak_anotation, y_pos + h * 0.15)

        x_boundary_min = x_boundary_min - w * 0.01
        x_boundary_max = x_boundary_max + w * 0.02
        if self.core.ncl == '13C':
            x_boundary_min = min(x_boundary_min, -10.0)
            x_boundary_max = max(x_boundary_max, 195.0)
        elif self.core.ncl == '1H':
            x_boundary_min = min(x_boundary_min, -0.2)
            x_boundary_max = max(x_boundary_max, 9.0)
        plt.xlim(
            x_boundary_max,
            x_boundary_min,
        )
        y_boundary_max = min(y_boundary_max, highest_peak_anotation)
        
        return y_boundary_max


    def __generate_info_box(self, plotlib):
        if not (self.core.is_sec or self.core.is_dsc):
            return
        core_dic = self.core.dic
        result = []
        if self.core.is_sec:
            sec_data_key = ['MN', 'MW', 'MP', 'D']
            for key in sec_data_key:
                dic_value = core_dic.get(key, [])
                key_str = f"{key}={dic_value[0]}" if len(dic_value) > 0 else None
                if key_str is not None:
                    result.append(key_str)
        else:
            dsc_meta_data = self.core.params.get('dsc_meta_data', None)
            melting_point, tg_value = '', ''
            if dsc_meta_data is not None:
                melting_point = dsc_meta_data.get('meltingPoint', '')
                tg_value = dsc_meta_data.get('tg', '')
            else:
                melting_point_arr = self.core.dic.get('MELTINGPOINT', [''])
                tg_value_arr = self.core.dic.get('TG', [''])
                melting_point = melting_point_arr[0]
                tg_value = tg_value_arr[0]
            
            melting_point_str = f"MELTING POINT={melting_point}"
            tg_str = f"TG={tg_value}"
            result.extend([melting_point_str, tg_str])

        info_str = '\n'.join(result)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax = plotlib.gca()
        ax.text(0.05, 0.95, info_str, fontsize=14, verticalalignment='top', bbox=props, transform=ax.transAxes)


    def __prepare_metadata_info_for_csv(self, csv_writer: csv.DictWriter):
        csv_writer.writerow({
            'Anodic E(V)': 'Measurement type',
            'Cathodic E(V)': 'Cyclic Voltammetry',
        })
        csv_writer.writerow({
            'Anodic E(V)': 'Measurement type ID',
        })
        csv_writer.writerow({
            'Anodic E(V)': 'Sample ID',
        })
        csv_writer.writerow({
            'Anodic E(V)': 'Analysis ID',
        })
        csv_writer.writerow({
            'Anodic E(V)': 'Dataset ID',
        })
        csv_writer.writerow({
            'Anodic E(V)': 'Dataset name',
        })
        csv_writer.writerow({
            'Anodic E(V)': 'Link to sample',
        })
        csv_writer.writerow({
        })

    def tf_csv(self):
        if self.core.is_cyclic_volta == False:
            return None
        tf_csv = tempfile.NamedTemporaryFile(suffix='.csv')
        tf_csv.flush()
        tf_csv.seek(0)

        listMaxMinPeaks = self.core.max_min_peaks_table
        if self.core.params['list_max_min_peaks'] is not None:
            listMaxMinPeaks = self.core.params['list_max_min_peaks']

        with open(tf_csv.name, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['Anodic E(V)', 'Anodic I(A)', 'Cathodic E(V)', 'Cathodic I(A)', 'I lambda0(A)', 'I ratio', 'E1/2(V)', 'Delta Ep(mV)']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            self.__prepare_metadata_info_for_csv(writer)

            writer.writeheader()

            x_max = np.max(self.core.xs)
            arr_max_x_indices = np.where(self.core.xs == x_max)
            y_pecker = 0
            if (arr_max_x_indices[0][0]):
                idx = arr_max_x_indices[0][0]
                y_pecker = self.core.ys[idx]

            for peak in listMaxMinPeaks:
                max_peak, min_peak, e12 = None, None, None
                if 'max' in peak:
                    max_peak = peak['max']
                if 'min' in peak:
                    min_peak = peak['min']
                if 'e12' in peak:
                    e12 = peak['e12']
                x_max_peak, y_max_peak = self.__get_xy_of_peak(max_peak)
                x_min_peak, y_min_peak = self.__get_xy_of_peak(min_peak)

                if (x_max_peak == '' and x_min_peak == ''):
                    # ignore if missing both max and min peak
                    continue

                if (x_max_peak == '' or x_min_peak == ''):
                    delta = ''
                else:
                    delta = abs(x_max_peak - x_min_peak) * 1000

                x_pecker = ''
                # calculate ratio
                if (y_min_peak == '' or y_max_peak == ''):
                    ratio = ''
                else:
                    if 'pecker' in peak and peak['pecker'] is not None:
                        pecker = peak['pecker']
                        x_pecker = pecker['x']
                        y_pecker = pecker['y']
                    first_expr = abs(y_min_peak) / abs(y_max_peak)
                    second_expr = 0.485 * abs(y_pecker) / abs(y_max_peak)
                    ratio = first_expr + second_expr + 0.086
                    if (y_pecker) == 0:
                        y_pecker = ''

                writer.writerow({
                    'Anodic E(V)': '{x_max}'.format(x_max=x_max_peak),
                    'Anodic I(A)': '{y_max}'.format(y_max=y_max_peak),
                    'Cathodic E(V)': '{x_min}'.format(x_min=x_min_peak),
                    'Cathodic I(A)': '{y_min}'.format(y_min=y_min_peak),
                    'I lambda0(A)': '{y_pecker}'.format(y_pecker=y_pecker),
                    'I ratio': '{ratio}'.format(ratio=ratio),
                    'E1/2(V)': '{e12}'.format(e12=e12),
                    'Delta Ep(mV)': '{delta}'.format(delta=delta)
                })
        return tf_csv

    def generate_nmrium(self, version=4, molfile_data=None):
        typ = self.core.typ
        if 'NMR' != typ:
            return None
        
        if version == 3:
            dic_data = self.__generate_nmrium_version_3(molfile_data=molfile_data)
        else:
            dic_data = self.__generate_nmrium_data(version=version, molfile_data=molfile_data)

        json_data = json.dumps(dic_data)
        
        tf_nmrium = tempfile.NamedTemporaryFile(suffix='.nmrium')
        tf_nmrium.write(bytes(json_data, 'UTF-8'))
        tf_nmrium.seek(0)
        return tf_nmrium
      
    def __generate_nmrium_version_3(self, molfile_data=None):
        dic_data = {'actionType': 'INITIATE', 'version': 3}
        spectra = self.__generate_nmrim_spectra()
        dic_data['spectra'] = spectra
        dic_data['molecules'] = self.__generate_molecules(molfile_data)
        return dic_data
      
    def __generate_nmrium_data(self, version, molfile_data=None):
        dic_data = {'version': version, 'data': {}, 'view': {}}
        spectra = self.__generate_nmrim_spectra()
        dic_data['data']['spectra'] = spectra
        dic_data['data']['molecules'] = self.__generate_molecules(molfile_data)
        return dic_data

    def __generate_nmrim_spectra(self):
        spectra = []
        spectra_id = str(uuid.uuid4())

        dic_spectra_data = {'id': spectra_id, 'source': {'jcampURL': None}}
        display_attr = {
            'name': spectra_id, 
            'color':'#C10020', 
            'isVisible':True, 
            'isPeaksMarkersVisible':True, 
            'isRealSpectrumVisible':True, 
            'isVisibleInDomain':True}
        dic_spectra_data['display'] = display_attr

        spectra_info = {'nucleus':self.core.ncl, 'isFid':False, 'dimension':1, 'isFt':True, 'type': 'NMR SPECTRUM'}
        dic_spectra_data['info'] = spectra_info
        dic_spectra_data['meta'] = self.core.dic

        x_values = np.flip(self.core.xs)
        y_values = np.flip(self.core.ys)
        dic_data_points = {'x':x_values.tolist(), 're':y_values.tolist(), 'im':None}
        dic_spectra_data['data'] = dic_data_points

        peaks = self.__generate_nmrim_peaks()
        dic_spectra_data['peaks'] = peaks

        integrals = self.__generate_nmrim_integrals()
        dic_spectra_data['integrals'] = integrals

        ranges = self.__generate_nmrim_ranges()
        dic_spectra_data['ranges'] = ranges

        spectra = [dic_spectra_data]

        return spectra

    def __generate_nmrim_peaks(self):
        dic_peaks = {'values':[], 'options':{}}

        x_peaks = []
        y_peaks = []
        if self.core.edit_peaks:
            x_peaks = self.core.edit_peaks['x']
            y_peaks = self.core.edit_peaks['y']
        elif self.core.auto_peaks:
            x_peaks = self.core.auto_peaks['x']
            y_peaks = self.core.auto_peaks['y']
        
        if len(x_peaks) != len(y_peaks):
            return dic_peaks

        for idx in range(len(x_peaks)):
            x = x_peaks[idx]
            y = y_peaks[idx]
            peak_id = str(uuid.uuid4())
            peak = {'id':peak_id, 'x':x, 'y':y, 'originalX': x}
            dic_peaks['values'].append(peak)

        return dic_peaks

    def __generate_nmrim_integrals(self):
        dic_integrals = {'values':[], 'options':{'isSumConstant':True, 'sumAuto':True, 'sum':100}}
        
        refShift, refArea = self.refShift, self.refArea
        if (len(self.all_itgs) == 0 and len(self.core.itg_table) > 0 and not self.core.params['integration'].get('edited') and ('originStack' not in self.core.params['integration'])):
            core_itg_table = self.core.itg_table[0]
            itg_table = core_itg_table.split('\n')
            for itg in itg_table:
                clear_itg = itg.replace('(', '')
                clear_itg = clear_itg.replace(')', '')
                split_itg = clear_itg.split(',')
                if (len(split_itg) > 2):
                    xLStr = split_itg[0].strip()
                    xUStr = split_itg[1].strip()
                    areaStr = split_itg[2].strip()
                    self.all_itgs.append({'xL': float(xLStr), 'xU': float(xUStr), 'area': float(areaStr)})    # noqa: E501
        for itg in self.all_itgs:
            xL, xU, area = itg['xL'] - refShift, itg['xU'] - refShift, itg['area'] * refArea    # noqa: E501
            iL, iU = get_curve_endpoint(self.core.xs, self.core.ys, xL, xU)
            cxs = self.core.xs[iL:iU]
            cys = self.core.ys[iL:iU]
            
            re_cxs = np.flip(cxs)
            re_cys = np.flip(cys)
            
            integral_id = str(uuid.uuid4())
            
            absolute_value = cal_xyIntegration(xs=re_cxs, ys=re_cys)

            integral = {'id':integral_id, 'originFrom':re_cxs[0], 'originTo':re_cxs[len(re_cxs)-1], 'from':re_cxs[0], 'to':re_cxs[len(re_cxs)-1], 'kind':'signal', 'absolute':absolute_value, 'integral':area*100}

            dic_integrals['values'].append(integral)
        return dic_integrals

    def __generate_nmrim_ranges(self):
        dic_ranges = {'values':[], 'options':{'isSumConstant':True, 'sumAuto':False, 'sum':100}}

        refShift, refArea = self.refShift, self.refArea
        if (len(self.mpys) == 0 and len(self.core.mpy_itg_table) > 0 and not self.core.params['integration'].get('edited') and ('originStack' not in self.core.params['integration'])):
            core_mpy_pks_table = self.core.mpy_pks_table[0]
            mpy_pks_table = core_mpy_pks_table.split('\n')
            tmp_dic_mpy_peaks = {}
            for peak in mpy_pks_table:
                clear_peak = peak.replace('(', '')
                clear_peak = clear_peak.replace(')', '')
                split_peak = clear_peak.split(',')
                idx_peakStr = split_peak[0].strip()
                xStr = split_peak[1].strip()
                yStr = split_peak[2].strip()
                if idx_peakStr not in tmp_dic_mpy_peaks:
                    tmp_dic_mpy_peaks[idx_peakStr] = []

                tmp_dic_mpy_peaks[idx_peakStr].append({'x': float(xStr), 'y': float(yStr)})

            core_mpy_itg_table = self.core.mpy_itg_table[0]
            mpy_itg_table = core_mpy_itg_table.split('\n')
            for mpy in mpy_itg_table:
                clear_mpy = mpy.replace('(', '')
                clear_mpy = clear_mpy.replace(')', '')
                split_mpy = clear_mpy.split(',')
                mpy_item = { 'mpyType': '', 'xExtent': {'xL': 0.0, 'xU': 0.0}, 'yExtent': {'yL': 0.0, 'yU': 0.0}, 'peaks': [], 'area': 1.0 }
                if (len(split_mpy) > 7):
                    idxStr = split_mpy[0].strip()
                    xLStr = split_mpy[1].strip()
                    xUStr = split_mpy[2].strip()
                    mpy_item['xExtent']['xL'] = float(xLStr) + refShift
                    mpy_item['xExtent']['xU'] = float(xUStr) + refShift
                    yLStr = split_mpy[3].strip()
                    yUStr = split_mpy[4].strip()
                    mpy_item['yExtent']['yL'] = float(yLStr) + refShift
                    mpy_item['yExtent']['yU'] = float(yUStr) + refShift
                    areaStr = split_mpy[5].strip()
                    mpy_item['area'] = float(areaStr)
                    typeStr = split_mpy[6].strip()
                    mpy_item['mpyType'] = typeStr
                    mpy_item['peaks'] = tmp_dic_mpy_peaks[idxStr]
                    self.mpys.append(mpy_item)
        
        total_integration_absolute = 0.0
        multiplicity_to_be_processed = self.mpys
        
        for mpy in multiplicity_to_be_processed:
            xL, xU, typ, peaks = mpy['xExtent']['xL'] - refShift, mpy['xExtent']['xU'] - refShift, mpy['mpyType'], mpy['peaks']    # noqa: E501
            iL, iU = get_curve_endpoint(self.core.xs, self.core.ys, xL, xU)
            cxs = self.core.xs[iL:iU]
            cys = self.core.ys[iL:iU]

            re_cxs = np.flip(cxs)
            re_cys = np.flip(cys)

            ranges_id = str(uuid.uuid4())
            mpy['ranges_id'] = ranges_id

            absolute_value = cal_xyIntegration(xs=re_cxs, ys=re_cys)
            mpy['absolute_value'] = absolute_value
            total_integration_absolute += absolute_value

            orgin_x_from, orgin_x_to = re_cxs[0], re_cxs[len(re_cxs)-1]
            mpy['orgin_x_from'] = orgin_x_from
            mpy['orgin_x_to'] = orgin_x_to
            
        for mpy in multiplicity_to_be_processed:
            typ, peaks = mpy['mpyType'], mpy['peaks']    # noqa: E501
            ranges_id, absolute_value = mpy['ranges_id'], mpy['absolute_value']
            orgin_x_from, orgin_x_to = mpy['orgin_x_from'], mpy['orgin_x_to']

            signal_id = str(uuid.uuid4())
            signal_delta = calc_mpy_center(mpy['peaks'], refShift, mpy['mpyType'])

            integration_value = (absolute_value / total_integration_absolute)*100

            signal_item = {'id':signal_id,'originDelta':signal_delta, 'delta':signal_delta, 'kind':'signal', 'integration':integration_value, 'multiplicity':typ, 'peaks':peaks}

            rang_item = {'id':ranges_id, 'originFrom':orgin_x_from, 'originTo':orgin_x_to, 'from':orgin_x_from, 'to':orgin_x_to, 'kind':'signal', 'absolute':absolute_value, 'integration':integration_value, 'signals':[signal_item]}
            dic_ranges['values'].append(rang_item)

        return dic_ranges
    
    def __generate_molecules(self, molfile_data):
        if molfile_data is None:
          return []
        
        molecule_id = str(uuid.uuid4())
        return [
            {
                'id': molecule_id,
                'label': 'P1',
                'molfile': molfile_data
            }
        ]
