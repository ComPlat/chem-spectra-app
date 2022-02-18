import matplotlib
matplotlib.use('Agg')

import tempfile  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.path as mpath  # noqa: E402
import re # noqa: E402
import numpy as np # noqa: E402

from chem_spectra.lib.composer.base import (  # noqa: E402
    extrac_dic, calc_npoints, BaseComposer
)
from chem_spectra.lib.shared.calc import (  # noqa: E402
    calc_mpy_center, calc_ks, get_curve_endpoint, cal_slope
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
        is_valid = ('13c' in nuc_modf) or ('1h' in nuc_modf) or ('19f' in nuc_modf) or ('31p' in nuc_modf) or ('15n' in nuc_modf) or ('29si' in nuc_modf)
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
        x_max, x_min = self.core.boundary['x']['max'], self.core.boundary['x']['min']
        xlim_left, xlim_right = [x_min, x_max] if (self.core.is_tga or self.core.is_uv_vis or self.core.is_hplc_uv_vis or self.core.is_xrd) else [x_max, x_min]
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
        if self.core.edit_peaks:
            plt.plot(
                self.core.edit_peaks['x'],
                self.core.edit_peaks['y'],
                'rd',
                marker=marker,
                markersize=50,
            )
        elif self.core.auto_peaks:
            plt.plot(
                self.core.auto_peaks['x'],
                self.core.auto_peaks['y'],
                'rd',
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
                  self.all_itgs.append({'xL': float(xLStr), 'xU': float(xUStr), 'area': float(areaStr)})
        for itg in self.all_itgs:
            # integration marker
            xL, xU, area = itg['xL'] - refShift, itg['xU'] - refShift, itg['area'] * refArea
            # integration curve
            ks = calc_ks(self.core.ys, y_max, h)
            iL, iU = get_curve_endpoint(self.core.xs, self.core.ys, xL, xU)
            ref = ks[iL]
            cxs = self.core.xs[iL:iU]

            if self.core.is_hplc_uv_vis:
                # fill area under curve
                cys = self.core.ys[iL:iU]
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
            else:
                # display integration
                plt.plot([xL, xU], [itg_h, itg_h], color='#228B22')
                plt.plot([xL, xL], [itg_h + h * 0.01, itg_h - h * 0.01], color='#228B22')
                plt.plot([xU, xU], [itg_h + h * 0.01, itg_h - h * 0.01], color='#228B22')
                plt.text((xL + xU) / 2, itg_h + h * 0.015, '{:0.2f}'.format(area), color='#228B22', ha='center', size=12)
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
            xL, xU, area, typ, peaks = mpy['xExtent']['xL'] - refShift, mpy['xExtent']['xU'] - refShift, mpy['area'] * refArea, mpy['mpyType'], mpy['peaks']
            plt.plot([xL, xU], [mpy_h, mpy_h], color='#DA70D6')
            plt.plot([xL, xL], [mpy_h + h * 0.01, mpy_h - h * 0.01], color='#DA70D6')
            plt.plot([xU, xU], [mpy_h + h * 0.01, mpy_h - h * 0.01], color='#DA70D6')
            plt.text((xL + xU) / 2, mpy_h - h * 0.1, '({})'.format(typ), color='#DA70D6', ha='center', size=12)
            plt.text((xL + xU) / 2, mpy_h - h * 0.06, '{:0.3f}'.format(calc_mpy_center(mpy['peaks'], refShift, mpy['mpyType'])), color='#DA70D6', ha='center', size=12)
            for p in peaks:
                x = p['x']
                plt.plot([x - refShift, x - refShift], [mpy_h, mpy_h + h * 0.05], color='#DA70D6')

        # PLOT label
        if (self.core.is_xrd):
            waveLength = self.core.params['waveLength']
            label = "X ({}), WL={} nm".format(self.core.label['x'], waveLength['value'], waveLength['unit'])
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
