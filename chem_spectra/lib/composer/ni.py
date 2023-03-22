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
import csv

from chem_spectra.lib.composer.base import (  # noqa: E402, F401
    extrac_dic, calc_npoints, BaseComposer
)
from chem_spectra.lib.shared.calc import (  # noqa: E402
    calc_mpy_center, calc_ks, get_curve_endpoint, cal_slope, cal_xyIntegration
)


TEXT_SPECTRUM_ORIG = '$$ === CHEMSPECTRA SPECTRUM ORIG ===\n'
TEXT_INTEGRATION = '$$ === CHEMSPECTRA INTEGRATION ===\n'
TEXT_MULTIPLICITY = '$$ === CHEMSPECTRA MULTIPLICITY ===\n'


class NIComposer(BaseComposer):
    def __init__(self, core):
        super().__init__(core)
        self.title = core.fname
        self.meta = self.__compose()

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
        return [
            '##.OBSERVE FREQUENCY={}\n'.format(
                extrac_dic(self.core, '.OBSERVEFREQUENCY')
            ),
            '##.OBSERVE NUCLEUS={}\n'.format(
                self.__get_nucleus()
            ),
            '##SPECTROMETER/DATA SYSTEM={}\n'.format(
                extrac_dic(self.core, 'SPECTROMETER/DATASYSTEM')
            ),
            '##.SHIFT REFERENCE={}\n'.format(
                extrac_dic(self.core, '.SHIFTREFERENCE')
            ),
            '##.SOLVENT NAME={}\n'.format(
                extrac_dic(self.core, '.SOLVENTNAME')
            ),
            '##.PULSE SEQUENCE={}\n'.format(
                extrac_dic(self.core, '.PULSESEQUENCE')
            ),
        ]

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

    def __get_xy_of_peak(self, peak):
        if peak is None:
            return '', ''
        x = peak['x'] or ''
        y = peak['y'] or ''
        return x, y

    def __gen_cyclic_voltammetry_data_peaks(self):
        content = ['##$CSCYCLICVOLTAMMETRYDATA=\n']
        if self.core.is_cyclic_volta:
            listMaxMinPeaks = self.core.max_min_peaks_table
            if self.core.params['list_max_min_peaks'] is not None:
                listMaxMinPeaks = self.core.params['list_max_min_peaks']

            x_max = np.max(self.core.xs)
            arr_max_x_indices = np.where(self.core.xs == x_max)
            y_pecker = 0
            if (arr_max_x_indices[0][0]):
                idx = arr_max_x_indices[0][0]
                y_pecker = self.core.ys[idx]

            for peak in listMaxMinPeaks:
                max_peak, min_peak = peak['max'], peak['min']
                x_max_peak, y_max_peak = self.__get_xy_of_peak(max_peak)
                x_min_peak, y_min_peak = self.__get_xy_of_peak(min_peak)

                if (x_max_peak == '' and x_min_peak == ''):
                    # ignore if missing both max and min peak
                    continue

                if (x_max_peak == '' or x_min_peak == ''):
                    delta = ''
                else:
                    delta = abs(x_max_peak - x_min_peak)

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

                content.append(
                    '({x_max}, {y_max}, {x_min}, {y_min}, {ratio}, {delta}, {x_pecker}, {y_pecker})\n'.format(x_max=x_max_peak, y_max=y_max_peak, x_min=x_min_peak, y_min=y_min_peak, ratio=ratio, delta=delta, x_pecker=x_pecker, y_pecker=y_pecker)  # noqa: E501
                )

        return content

    def __compose(self):
        meta = []
        meta.extend(self.gen_headers_root())

        meta.extend(self.__gen_headers_spectrum_orig())
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
            meta.extend(self.__gen_cyclic_voltammetry_data_peaks())
        meta.extend(self.gen_ending())

        meta.extend(self.gen_headers_peaktable_edit())
        meta.extend(self.gen_edit_peaktable())
        meta.extend(self.gen_ending())

        meta.extend(self.gen_headers_peaktable_auto())
        meta.extend(self.gen_auto_peaktable())
        meta.extend(self.gen_ending())

        meta.extend(self.gen_ending())
        return meta

    def __plt_nbins(self):
        typ = self.core.typ
        if 'NMR' == typ:
            return 20
        elif 'INFRARED' == typ:
            return 20
        elif 'RAMAN' == typ:
            return 20
        elif 'UVVIS' == typ:
            return 20
        elif 'THERMOGRAVIMETRIC ANALYSIS' == typ:
            return 20
        elif 'MS' == typ:
            return 20
        return 20

    def __fakto(self):
        typ = self.core.typ
        if 'NMR' == typ:
            return 1
        elif 'MS' == typ:
            return 1
        elif 'INFRARED' == typ:
            return -1
        elif 'RAMAN' == typ:
            return 1
        elif 'UVVIS' == typ:
            return 1
        elif 'THERMOGRAVIMETRIC ANALYSIS' == typ:
            return 1
        return 1

    def tf_img(self):
        plt.rcParams['figure.figsize'] = [16, 9]
        plt.rcParams['font.size'] = 14

        # PLOT data
        plt.plot(self.core.xs, self.core.ys)
        x_max, x_min = self.core.boundary['x']['max'], self.core.boundary['x']['min']   # noqa: E501
        xlim_left, xlim_right = [x_min, x_max] if (self.core.is_tga or self.core.is_uv_vis or self.core.is_hplc_uv_vis or self.core.is_xrd or self.core.is_cyclic_volta or self.core.is_sec) else [x_max, x_min]    # noqa: E501
        plt.xlim(xlim_left, xlim_right)
        y_max, y_min = np.max(self.core.ys), np.min(self.core.ys)
        h = y_max - y_min
        plt.ylim(
            y_min - h * 0.2,
            y_max + h * 0.2,
        )
        # PLOT peaks
        faktor = self.__fakto()
        path_data = [
            (mpath.Path.MOVETO, (0, faktor * 5)),
            (mpath.Path.LINETO, (0, faktor * 20)),
        ]
        codes, verts = zip(*path_data)
        marker = mpath.Path(verts, codes)
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
        if self.core.is_cyclic_volta:
            x_peaks = []
            y_peaks = []
            formatter = ticker.ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((-1, 1))
            plt.gca().yaxis.set_major_formatter(formatter)

            listMaxMinPeaks = []
            if self.core.params['list_max_min_peaks'] is not None:
                listMaxMinPeaks = self.core.params['list_max_min_peaks']

            for peak in listMaxMinPeaks:
                max_peak, min_peak = peak['max'], peak['min']
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

        plt.plot(
            x_peaks,
            y_peaks,
            'rd',
            marker=marker,
            markersize=50,
        )

        plt.plot(
            x_peckers,
            y_peckers,
            'gd',
            marker=marker,
            markersize=50,
        )

        # ----- PLOT integration -----
        refShift, refArea = self.refShift, self.refArea
        itg_h = y_max + h * 0.1
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
                plt.plot([xL, xU], [itg_h, itg_h], color='#228B22')
                plt.plot([xL, xL], [itg_h + h * 0.01, itg_h - h * 0.01], color='#228B22')   # noqa: E501
                plt.plot([xU, xU], [itg_h + h * 0.01, itg_h - h * 0.01], color='#228B22')   # noqa: E501
                plt.text((xL + xU) / 2, itg_h + h * 0.015, '{:0.2f}'.format(area), color='#228B22', ha='center', size=12)   # noqa: E501
                cys = (ks[iL:iU] - ref) * 1.5 + (y_max - h * 0.4)
                # if self.core.is_uv_vis:
                #     cys = (ref - ks[iL:iU]) * 0.5 + (y_max - h * 0.4)
                plt.plot(cxs, cys, color='#228B22')

        # ----- PLOT multiplicity -----
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
        mpy_h = y_min - h * 0.08
        for mpy in self.mpys:
            xL, xU, area, typ, peaks = mpy['xExtent']['xL'] - refShift, mpy['xExtent']['xU'] - refShift, mpy['area'] * refArea, mpy['mpyType'], mpy['peaks']    # noqa: E501
            plt.plot([xL, xU], [mpy_h, mpy_h], color='#DA70D6')
            plt.plot([xL, xL], [mpy_h + h * 0.01, mpy_h - h * 0.01], color='#DA70D6')   # noqa: E501
            plt.plot([xU, xU], [mpy_h + h * 0.01, mpy_h - h * 0.01], color='#DA70D6')   # noqa: E501
            plt.text((xL + xU) / 2, mpy_h - h * 0.1, '({})'.format(typ), color='#DA70D6', ha='center', size=12)  # noqa: E501
            plt.text((xL + xU) / 2, mpy_h - h * 0.06, '{:0.3f}'.format(calc_mpy_center(mpy['peaks'], refShift, mpy['mpyType'])), color='#DA70D6', ha='center', size=12)  # noqa: E501
            for p in peaks:
                x = p['x']
                plt.plot([x - refShift, x - refShift], [mpy_h, mpy_h + h * 0.05], color='#DA70D6')  # noqa: E501

        # PLOT label
        if (self.core.is_xrd):
            waveLength = self.core.params['waveLength']
            label = "X ({}), WL={} nm".format(self.core.label['x'], waveLength['value'], waveLength['unit'])    # noqa: E501
            plt.xlabel((label), fontsize=18)
        else:
            plt.xlabel("X ({})".format(self.core.label['x']), fontsize=18)
        plt.ylabel("Y ({})".format(self.core.label['y']), fontsize=18)
        plt.locator_params(nbins=self.__plt_nbins())
        plt.grid(False)
        # Save
        tf_img = tempfile.NamedTemporaryFile(suffix='.png')
        plt.savefig(tf_img, format='png')
        tf_img.seek(0)
        plt.clf()
        plt.cla()
        return tf_img

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
            # fieldnames = ['Max', 'Min', 'I λ0', 'I ratio', 'Pecker']
            fieldnames = ['Max x', 'Max y', 'Min x', 'Min y', 'Delta Ep', 'I lambda0', 'I ratio']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            x_max = np.max(self.core.xs)
            arr_max_x_indices = np.where(self.core.xs == x_max)
            y_pecker = 0
            if (arr_max_x_indices[0][0]):
                idx = arr_max_x_indices[0][0]
                y_pecker = self.core.ys[idx]

            for peak in listMaxMinPeaks:
                max_peak, min_peak = peak['max'], peak['min']
                x_max_peak, y_max_peak = self.__get_xy_of_peak(max_peak)
                x_min_peak, y_min_peak = self.__get_xy_of_peak(min_peak)

                if (x_max_peak == '' and x_min_peak == ''):
                    # ignore if missing both max and min peak
                    continue

                if (x_max_peak == '' or x_min_peak == ''):
                    delta = ''
                else:
                    delta = abs(x_max_peak - x_min_peak)

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
                    'Max x': '{x_max}'.format(x_max=x_max_peak),
                    'Max y': '{y_max}'.format(y_max=y_max_peak),
                    'Min x': '{x_min}'.format(x_min=x_min_peak),
                    'Min y': '{y_min}'.format(y_min=y_min_peak),
                    'Delta Ep': '{y_pecker}'.format(y_pecker=y_pecker),
                    'I lambda0': '{ratio}'.format(ratio=ratio),
                    'I ratio': '{delta}'.format(delta=delta)
                })
        return tf_csv

    def generate_nmrium(self):
        typ = self.core.typ
        if 'NMR' != typ:
            return None
        
        dic_data = {'actionType': 'INITIATE', 'version': 3}
        spectra = self.__generate_nmrim_spectra()
        dic_data['spectra'] = spectra

        json_data = json.dumps(dic_data)
        
        tf_nmrium = tempfile.NamedTemporaryFile(suffix='.nmrium')
        tf_nmrium.write(bytes(json_data, 'UTF-8'))
        tf_nmrium.seek(0)
        return tf_nmrium

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
            peak = {'id':peak_id, 'x':x, 'y':y}
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
