import numpy as np
from scipy import signal

from chem_spectra.lib.converter.datatable import DatatableModel


THRESHOLD_IR = 0.93
THRESHOLD_RAMAN = 0.07
THRESHOLD_NMR = 0.005
THRESHOLD_MS = 0.05


class JcampNIConverter:  # nmr & IR
    def __init__(self, base):
        self.params = base.params
        self.dic = base.dic
        self.data = base.data
        self.datatypes = base.datatypes
        self.datatype = base.datatype
        self.title = base.title
        self.typ = base.typ
        self.is_em_wave = base.is_em_wave
        self.is_ir = base.is_ir
        # - - - - - - - - - - -
        self.fname = base.fname
        self.target_idx = self.__index_target()
        self.block_count = self.__count_block()
        self.threshold = self.__thres()
        self.ncl = self.__ncl()
        self.obs_freq = self.__set_obs_freq()
        self.factor = self.__set_factor()
        self.x_unit = self.__set_x_unit()
        self.ys = self.__read_ys()
        self.xs = self.__read_xs()
        self.boundary = self.__find_boundary()
        self.label = self.__set_label()
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
        return 0.5

    def __ncl(self):
        try:
            ncls = self.dic['.OBSERVENUCLEUS']
            if '^1H' in ncls:
                return '1H'
            elif '^13C' in ncls:
                return '13C'
            elif '^19F' in ncls:
                return '19F'
        except: # noqa
            pass
        return ''

    def __index_target(self):
        target_topics = [
            'NMR SPECTRUM', 'NMRSPECTRUM',
            'INFRARED SPECTRUM', 'RAMAN SPECTRUM',
            'MASS SPECTRUM'
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
                beg_pt = float(self.dic['FIRSTX'][idx])
                end_pt = float(self.dic['LASTX'][idx])
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
                'x': self.dic['XFACTOR'][0],
                'y': self.dic['YFACTOR'][0],
            }
        except:  # noqa
            pass

        return factor

    def __set_x_unit(self):
        x_unit = None

        try:
            x_unit = self.dic['XUNITS'][0].upper()
        except:  # noqa
            pass

        return x_unit

    def __set_auto_peaks(self, peak_idxs):
        auto_x = []
        auto_y = []
        for idx in peak_idxs:
            auto_x.append(self.xs[idx])
            auto_y.append(self.ys[idx])
        if len(auto_x) == 0:
            return
        self.auto_peaks = {'x': auto_x, 'y': auto_y}

    def __read_auto_peaks(self):
        if self.params['clear']:
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
                if len(self.dic['PEAKTABLE']):
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
        if self.params['clear']:
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

    def __exec_peak_picking_logic(self):
        max_y = np.max(self.ys)
        height = self.threshold * max_y

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
        within_limit = False
        while not within_limit:
            peak_idxs = self.__exec_peak_picking_logic()
            if peak_idxs.shape[0] <= 100:
                within_limit = True
                self.__set_auto_peaks(peak_idxs)
            else:
                if self.is_ir:
                    self.threshold *= 0.9
                else:
                    self.threshold *= 1.5

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
