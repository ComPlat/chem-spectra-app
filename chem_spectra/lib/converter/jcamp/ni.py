import numpy as np
from scipy import signal

from chem_spectra.lib.converter.datatable import DatatableModel
from chem_spectra.lib.shared.calc import to_float
from chem_spectra.lib.shared.equal_space import make_equal_space_1d


THRESHOLD_IR = 0.93
THRESHOLD_RAMAN = 0.07
THRESHOLD_NMR = 0.005
THRESHOLD_MS = 0.05
THRESHOLD_UVVIS = 0.05
THRESHOLD_TGA = 1.05

class JcampNIConverter:  # nmr & IR
    def __init__(self, base):
        self.params = base.params
        self.dic = base.dic
        self.data = make_equal_space_1d(base)
        self.datatypes = base.datatypes
        self.datatype = base.datatype
        self.title = base.title
        self.typ = base.typ
        self.is_em_wave = base.is_em_wave
        self.is_ir = base.is_ir
        self.non_nmr = base.non_nmr
        self.ncl = base.ncl
        self.solv_peaks = base.solv_peaks
        # - - - - - - - - - - -
        self.fname = base.fname
        self.target_idx = self.__index_target()
        self.block_count = self.__count_block()
        self.threshold = self.__thres()
        self.obs_freq = self.__set_obs_freq()
        self.factor = self.__set_factor()
        self.x_unit = self.__set_x_unit()
        self.ys = self.__read_ys()
        self.xs = self.__read_xs()
        self.__set_first_last_xs()
        self.clear = self.__refresh_solvent()
        self.boundary = self.__find_boundary()
        self.label = self.__set_label()
        self.simu_peaks = base.simu_peaks
        self.auto_peaks = None
        self.edit_peaks = None
        self.itg_table = []
        self.mpy_itg_table = []
        self.mpy_pks_table = []
        self.datatable = self.__set_datatable()
        self.__read_peak_from_file()
        self.__read_integration_from_file()
        self.__read_multiplicity_from_file()

    def __thres(self):
        dt = self.datatype
        if 'NMR SPECTRUM' == dt:
            return THRESHOLD_NMR
        elif 'NMRSPECTRUM' == dt:  # MNova
            return THRESHOLD_NMR
        elif 'INFRARED SPECTRUM' == dt:
            return THRESHOLD_IR
        elif 'RAMAN SPECTRUM' == dt:
            return THRESHOLD_RAMAN
        elif 'MASS SPECTRUM' == dt:
            return THRESHOLD_MS
        elif 'UV/VIS SPECTRUM' == dt:
            return THRESHOLD_UVVIS
        elif 'UV-VIS' == dt:
            return THRESHOLD_UVVIS
        elif 'THERMOGRAVIMETRIC ANALYSIS' == dt:
            return THRESHOLD_TGA
        return 0.5

    def __index_target(self):
        target_topics = [
            'NMR SPECTRUM', 'NMRSPECTRUM',
            'INFRARED SPECTRUM', 'RAMAN SPECTRUM',
            'MASS SPECTRUM', 'UV/VIS SPECTRUM', 'UV-VIS',
            'THERMOGRAVIMETRIC ANALYSIS'
        ]
        for tp in target_topics:
            if tp in self.datatypes:
                idx = self.datatypes.index(tp)

        if 'LINK' in self.datatypes:
            idx -= 1

        return idx

    def __count_block(self):
        count = 1
        try:
            count = int(self.dic['BLOCKS'][0])
        except:  # noqa
            pass

        return count

    def __read_xs(self):  # TBD
        beg_pt = None
        end_pt = None
        idx = self.target_idx

        if beg_pt is None:
            try:
                obs_freq = self.obs_freq
                shift = float(self.dic['$OFFSET'][idx])
                beg_pt = float(
                    self.dic['FIRST'][idx].replace(' ', '').split(',')[0]
                ) / obs_freq
                end_pt = float(
                    self.dic['LAST'][idx].replace(' ', '').split(',')[0]
                ) / obs_freq
                shift = beg_pt - shift
                beg_pt = beg_pt - shift
                end_pt = end_pt - shift
            except:  # noqa
                pass

        if beg_pt is None:  # MNova
            try:
                obs_freq = self.obs_freq
                beg_pt = float(
                    self.dic['FIRST'][idx].replace(' ', '').split(',')[0]
                ) / obs_freq
                end_pt = float(
                    self.dic['LAST'][idx].replace(' ', '').split(',')[0]
                ) / obs_freq
            except:  # noqa
                pass

        if beg_pt is None:
            try:
                beg_pt = to_float(self.dic['FIRSTX'][idx])
                end_pt = to_float(self.dic['LASTX'][idx])
            except:  # noqa
                pass

        if self.is_em_wave and beg_pt < end_pt:
            buf = beg_pt
            beg_pt = end_pt
            end_pt = buf
            self.ys = self.ys[::-1]

        num_pt = self.ys.shape[0]
        x = np.linspace(
            beg_pt + self.params['delta'],
            end_pt + self.params['delta'],
            num=num_pt,
            endpoint=True
        )

        if self.x_unit == 'HZ':
            x = x / self.obs_freq

        return x

    def __read_ys(self):
        y = None
        try:
            y = self.data['real'][self.target_idx]
        except:  # noqa
            pass

        if y is None:
            try:
                y = self.data[:]
            except:  # noqa
                pass
        # transmission only # IR ABS vs TRANS
        if self.is_ir:
            y_median = np.median(y)
            y_max = np.max(y)
            if y_median < 0.5 * y_max:
                y = y_max - y

        return y

    def __find_boundary(self):
        return {
            'x': {
                'max': self.xs.max(),
                'min': self.xs.min(),
            },
            'y': {
                'max': self.ys.max(),
                'min': self.ys.min(),
            },
        }

    def __set_label(self):
        target = {'x': 'PPM', 'y': 'ARBITRARY'}
        try:
            x = self.dic['XUNITS'][self.target_idx]
            y = self.dic['YUNITS'][self.target_idx]
            x = 'PPM' if x.upper() == 'HZ' else x
            y = 'ARBITRARY' if y.upper() == 'ARBITRARYUNITS' else y
            target = {'x': x, 'y': y}
        except:  # noqa
            pass

        try:
            x, y, _ = self.dic['UNITS'][1].replace(' ', '').split(',')
            x = 'PPM' if x.upper() == 'HZ' else x
            y = 'ARBITRARY' if y.upper() == 'ARBITRARYUNITS' else y
            target = {'x': x, 'y': y}
        except:  # noqa
            pass

        if 'absorb' in target['y'].lower():  # IR ABS vs TRANS
            target['y'] = 'TRANSMITTANCE'

        return target

    def __set_obs_freq(self):
        obs_freq = None
        try:
            obs_freq = float(self.dic['.OBSERVEFREQUENCY'][0])
        except:  # noqa
            pass

        try:
            if obs_freq is None:
                obs_freq = float(self.dic['$SFO1'][0])
        except:  # noqa
            pass

        return obs_freq

    def __set_factor(self):
        factor = {
            'x': 1.0,
            'y': 1.0,
        }

        try:
            factor = {
                'x': to_float(self.dic['XFACTOR'][0]),
                'y': to_float(self.dic['YFACTOR'][0]),
            }
        except:  # noqa
            pass

        if factor['y'] == 1.0 and not isinstance(self.data, dict):
            factor['y'] = self.data.max() / 1000000.0

        return factor

    def __set_x_unit(self):
        x_unit = None

        try:
            x_unit = self.dic['XUNITS'][0].upper()
        except:  # noqa
            pass

        return x_unit

    def __read_auto_peaks(self):
        if self.params['clear'] or self.clear:
            return

        try:  # legacy
            auto_x = []
            auto_y = []
            pas = self.dic['PEAKASSIGNMENTS'][0].split('\n')[1:]
            for pa in pas:
                info = pa.replace('(', '').replace(')', '') \
                            .replace(' ', '').split(',')
                auto_x.append(float(info[0]))
                auto_y.append(float(info[1]))
            if len(auto_x) == 0:
                return
            self.auto_peaks = {'x': auto_x, 'y': auto_y}
        except:  # noqa
            pass

        try:  # mnova
            if self.auto_peaks is None:
                auto_x = []
                auto_y = []
                if len(self.dic['PEAKTABLE']) == 0:
                    return
                pas = self.dic['PEAKTABLE'][1].split('\n')[1:]
                for pa in pas:
                    info = pa.replace(' ', '').split(',')
                    auto_x.append(float(info[0]))
                    auto_y.append(float(info[1]))
                if len(auto_x) == 0:
                    return
                self.auto_peaks = {'x': auto_x, 'y': auto_y}
        except:  # noqa
            pass

    def __read_edit_peaks(self):
        if self.params['clear'] or self.clear:
            return

        try:  # legacy
            edit_x = []
            edit_y = []
            pas = self.dic['PEAKASSIGNMENTS'][1].split('\n')[1:]
            for pa in pas:
                info = pa.replace('(', '').replace(')', '') \
                            .replace(' ', '').split(',')
                edit_x.append(float(info[0]))
                edit_y.append(float(info[1]))
            if len(edit_x) == 0:
                return
            self.edit_peaks = {'x': edit_x, 'y': edit_y}
        except:  # noqa
            pass

        try:  # mnova
            if self.edit_peaks is None:
                edit_x = []
                edit_y = []
                pas = self.dic['PEAKTABLE'][0].split('\n')[1:]
                for pa in pas:
                    info = pa.replace(' ', '').split(',')
                    edit_x.append(float(info[0]))
                    edit_y.append(float(info[1]))
                if len(edit_x) == 0:
                    return
                self.edit_peaks = {'x': edit_x, 'y': edit_y}
        except:  # noqa
            pass

    def __parse_edit(self):
        edit_x = []
        edit_y = []
        for p in self.params['peaks_str'].split('#'):
            info = p.split(',')
            edit_x.append(float(info[0]))
            edit_y.append(float(info[1]))
        if len(edit_x) == 0:
            return
        self.edit_peaks = {'x': edit_x, 'y': edit_y}

    def __exec_peak_picking_logic(self, refresh_solvent=False):
        max_y = np.max(self.ys)
        height = 0.2 * max_y if refresh_solvent else self.threshold * max_y

        corr_data_ys = self.ys
        corr_height = height
        if self.is_ir:
            corr_data_ys = 1 - self.ys
            corr_height = 1 - height

        peak_idxs = signal.find_peaks(corr_data_ys, height=corr_height)[0]

        min_y = np.min(self.ys)
        if not self.is_ir and (max_y * 0.4 < -min_y):
            dept_corr_data_ys = 1 - self.ys
            dept_corr_height = height
            dept_peak_idxs = signal.find_peaks(dept_corr_data_ys, height=dept_corr_height)[0]
            peak_idxs = np.unique(np.concatenate((peak_idxs, dept_peak_idxs)))
        return peak_idxs

    def __run_auto_pick_peak(self):
        peak_idxs = self.__exec_peak_picking_logic()
        auto_peaks = [{'x': self.xs[idx], 'y': self.ys[idx]} for idx in peak_idxs]
        auto_peaks.sort(key=lambda d: d['y'], reverse=True)

        if self.is_ir:
            auto_peaks = auto_peaks[-100:]
        elif self.ncl == '13C':
            simu_length = len(self.simu_peaks)
            simu_length = simu_length if simu_length > 1 else 50
            auto_peaks = auto_peaks[:200]
            # rm solvent peaks
            edit_non_solv_peaks = []
            for peak in auto_peaks:
                not_solvent = True
                for u, v in self.solv_peaks:
                    if u < peak['x'] < v:
                        not_solvent = False
                if not_solvent:
                    edit_non_solv_peaks.append(peak)
            # as - 26.90 (range: 26.80 - 26.95)
            # bs, cs - 207.1 (207.0-207.2) + 30.9 (30.8 - 31.0)
            # ds, es, fs, gs - 60.4 (60.3 - 60.5) + 14.2 (14.1 - 14.3) + 21.1 (20.9 - 21.2) + 171.3 (171.2-171.4)
            # rm impurity peaks
            imp_as, imp_bs, imp_cs, imp_ds, imp_es, imp_fs, imp_gs, edit_pure_peaks = [], [], [], [], [], [], [], []
            i, capacity, l = 0, simu_length + 10, len(edit_non_solv_peaks)
            while i < capacity and i < l:
                target = edit_non_solv_peaks[i]
                if 26.80 <= target['x'] <= 27.0:
                    imp_as.append(target)
                elif 207.0 <= target['x'] <= 207.2:
                    imp_bs.append(target)
                elif 30.8 <= target['x'] <= 31.0:
                    imp_cs.append(target)
                elif 60.3 <= target['x'] <= 60.5:
                    imp_ds.append(target)
                elif 14.1 <= target['x'] <= 14.3:
                    imp_es.append(target)
                elif 20.9 <= target['x'] <= 21.2:
                    imp_fs.append(target)
                elif 171.1 <= target['x'] <= 171.4:
                    imp_gs.append(target)
                else:
                    edit_pure_peaks.append(target)
                    i += 1
                    continue
                i += 1
                capacity += 1
            if not (imp_bs and imp_cs):
                edit_pure_peaks = edit_pure_peaks + imp_bs + imp_cs
            if not (imp_ds and imp_es and imp_fs and imp_gs):
                edit_pure_peaks = edit_pure_peaks + imp_ds + imp_es + imp_fs + imp_gs
            edit_pure_peaks.sort(key=lambda d: d['y'], reverse=True)
            edit_peaks = edit_pure_peaks[:simu_length]
            # assign to edit_peaks
            edit_x = [peak['x'] for peak in edit_peaks]
            edit_y = [peak['y'] for peak in edit_peaks]
            self.edit_peaks = {'x': edit_x, 'y': edit_y}
            auto_peaks = auto_peaks[:100]
        elif self.ncl == '1H':
            auto_peaks = auto_peaks[:100]
            edit_non_solv_peaks = []
            for peak in auto_peaks:
                not_solvent = True
                for u, v in self.solv_peaks:
                    if u < peak['x'] < v:
                        not_solvent = False
                if not_solvent:
                    edit_non_solv_peaks.append(peak)
            edit_x = [peak['x'] for peak in edit_non_solv_peaks]
            edit_y = [peak['y'] for peak in edit_non_solv_peaks]
            self.edit_peaks = {'x': edit_x, 'y': edit_y}
        else:
            auto_peaks = auto_peaks[:100]

        auto_x = [peak['x'] for peak in auto_peaks]
        auto_y = [peak['y'] for peak in auto_peaks]
        self.auto_peaks = {'x': auto_x, 'y': auto_y}

    def __set_datatable(self):
        y_factor = self.factor and self.factor['y']
        y_factor = y_factor or 1.0
        return DatatableModel().encode(
            self.ys,
            y_factor
        )

    def __read_peak_from_file(self):
        self.__read_auto_peaks()
        self.__read_edit_peaks()
        if not self.auto_peaks or not self.params['delta'] == 0.0:
            self.__run_auto_pick_peak()
        if self.params['peaks_str']:
            self.__parse_edit()

    def __read_integration_from_file(self):
        target = self.dic.get('$OBSERVEDINTEGRALS')
        if target:
            self.itg_table = ['\n'.join(target[0].split('\n')[1:]), '\n']

    def __read_multiplicity_from_file(self):
        target1 = self.dic.get('$OBSERVEDMULTIPLETS')
        if target1:
            self.mpy_itg_table = target1
            self.mpy_itg_table.append('\n')
        target2 = self.dic.get('$OBSERVEDMULTIPLETSPEAKS')
        if target2:
            self.mpy_pks_table = target2
            self.mpy_pks_table.append('\n')

    def __refresh_solvent(self):
        if self.ncl == '13C':
            ref_name = (
                self.params['ref_name'] or
                self.dic.get('$CSSOLVENTNAME', [''])[0]
            )
            if ref_name and ref_name != '- - -':
                return
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            peak_idxs = self.__exec_peak_picking_logic(refresh_solvent=True)[:100]
            auto_peaks = [{'x': self.xs[idx], 'y': self.ys[idx], 'idx': idx} for idx in peak_idxs]
            auto_peaks.sort(key=lambda d: d['y'], reverse=True)
            # is acetone
            left, right = auto_peaks[0], auto_peaks[1]
            if left['x'] < right['x']: left, right = right, left
            diff = abs(left['x'] - right['x'])
            if 175 < diff < 177:
                self.dic['$CSSOLVENTNAME'] = ['Acetone-d6 (sep)']
                self.dic['$CSSOLVENTVALUE'] = ['29.920']
                self.dic['$CSSOLVENTX'] = ['0']
                self.solv_peaks = [(27.0, 33.0), (203.7, 209.7)]
                shift = 29.920 - right['x']
                self.xs = self.xs + shift
                return True # self.clear
            # is chloroform
            for hpk in auto_peaks[:10]:
                x_c = hpk['x']
                peaks = [p for p in auto_peaks if x_c - 2.0 < p['x'] < x_c + 2.0]
                if len(peaks) == 3:
                    pxs = sorted(map(lambda p: p['x'], peaks))
                    diff_one = abs(pxs[0] - pxs[1])
                    diff_two = abs(pxs[1] - pxs[2])
                    if 0.2 < diff_one < 0.6 and 0.2 < diff_two < 0.6:
                        self.dic['$CSSOLVENTNAME'] = ['Chloroform-d (t)']
                        self.dic['$CSSOLVENTVALUE'] = ['77.00']
                        self.dic['$CSSOLVENTX'] = ['0']
                        self.solv_peaks = [(74.0, 80.0)]
                        shift = 77.00 - pxs[1]
                        self.xs = self.xs + shift
                        return True # self.clear

        return False # self.clear

    def __set_first_last_xs(self):
        self.first_x = self.xs[0]
        self.last_x = self.xs[-1]