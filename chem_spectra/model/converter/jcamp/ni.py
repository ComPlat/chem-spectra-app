import nmrglue as ng
import numpy as np
from scipy import signal

from chem_spectra.model.converter.encoder import encode_datatable


THRESHOLD_IR = 0.93
THRESHOLD_NMR = 0.005
THRESHOLD_MS = 0.05


class JcampNIConverter(): # nmr & IR
    def __init__(self, base):
        self.params = base.params
        self.dic = base.dic
        self.data = base.data
        self.datatypes = base.datatypes
        self.datatype = base.datatype
        self.title = base.title
        self.typ = base.typ
        # - - - - - - - - - - -
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
        self.datatable = self.__set_datatable()
        self.__read_peak_from_file()


    def __thres(self):
        dt = self.datatype
        if 'NMR SPECTRUM' == dt:
            return THRESHOLD_NMR
        elif 'NMRSPECTRUM' == dt: # MNova
            return THRESHOLD_NMR
        elif 'INFRARED SPECTRUM' == dt:
            return THRESHOLD_IR
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
        except:
            pass
        return ''


    def __index_target(self):
        target_topics = ['NMR SPECTRUM', 'NMRSPECTRUM', 'INFRARED SPECTRUM', 'MASS SPECTRUM']
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
        except:
            pass

        return count


    def __read_xs(self): # TBD
        beg_pt = None
        end_pt = None
        idx = self.target_idx

        if beg_pt is None:
            try:
                obs_freq = self.obs_freq
                shift = float(self.dic['$OFFSET'][idx])
                beg_pt = float(self.dic['FIRST'][idx].replace(' ', '').split(',')[0]) / obs_freq
                end_pt = float(self.dic['LAST'][idx].replace(' ', '').split(',')[0]) / obs_freq
                shift = beg_pt - shift
                beg_pt = beg_pt - shift
                end_pt = end_pt - shift
            except:
                pass

        if beg_pt is None: # MNova
            try:
                obs_freq = self.obs_freq
                beg_pt = float(self.dic['FIRST'][idx].replace(' ', '').split(',')[0]) / obs_freq
                end_pt = float(self.dic['LAST'][idx].replace(' ', '').split(',')[0]) / obs_freq
            except:
                pass

        if beg_pt is None:
            try:
                beg_pt = float(self.dic['FIRSTX'][idx])
                end_pt = float(self.dic['LASTX'][idx])
            except:
                pass

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
        except:
            pass

        if y is None:
            try:
                y = self.data[:]
            except:
                pass

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
        try:
            x = self.dic['XUNITS'][self.target_idx]
            y = self.dic['YUNITS'][self.target_idx]
            x = 'PPM' if x.upper() == 'HZ' else x
            y = 'ARBITRARY' if y.upper() == 'ARBITRARYUNITS' else y
            return { 'x': x, 'y': y }
        except:
            pass

        try:
            x, y, _ = dic['UNITS'][1].replace(' ', '').split(',')
            x = 'PPM' if x.upper() == 'HZ' else x
            y = 'ARBITRARY' if y.upper() == 'ARBITRARYUNITS' else y
            return { 'x': x, 'y': y }
        except:
            pass

        return { 'x': 'PPM', 'y': 'ARBITRARY' }


    def __set_obs_freq(self):
        obs_freq = None
        try:
            obs_freq = float(self.dic['.OBSERVEFREQUENCY'][0])
        except:
            pass

        try:
            if obs_freq is None:
                obs_freq = float(self.dic['$SFO1'][0])
        except:
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
        except:
            pass

        return factor


    def __set_x_unit(self):
        x_unit = None

        try:
            x_unit = self.dic['XUNITS'][0].upper()
        except:
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
        self.auto_peaks = { 'x': auto_x, 'y': auto_y }


    def __read_auto_peaks(self):
        try: # legacy
            auto_x = []
            auto_y = []
            pas = self.dic['PEAKASSIGNMENTS'][0].split('\n')[1:]
            for pa in pas:
                info = pa.replace('(', '').replace(')', '').replace(' ', '').split(',')
                auto_x.append(float(info[0]))
                auto_y.append(float(info[1]))
            if len(auto_x) == 0:
                return
            self.auto_peaks = { 'x': auto_x, 'y': auto_y }
        except:
            pass

        try: # mnova
            if self.auto_peaks is None:
                auto_x = []
                auto_y = []
                if len(self.dic['PEAKTABLE']):
                    return
                pas = self.dic['PEAKTABLE'][-1].split('\n')[1:]
                for pa in pas:
                    info = pa.replace(' ', '').split(',')
                    auto_x.append(float(info[0]))
                    auto_y.append(float(info[1]))
                if len(auto_x) == 0:
                    return
                self.auto_peaks = { 'x': auto_x, 'y': auto_y }
        except:
            pass


    def __read_edit_peaks(self):
        try: # legacy
            edit_x = []
            edit_y = []
            pas = self.dic['PEAKASSIGNMENTS'][-1].split('\n')[1:]
            for pa in pas:
                info = pa.replace('(', '').replace(')', '').replace(' ', '').split(',')
                edit_x.append(float(info[0]))
                edit_y.append(float(info[1]))
            if len(edit_x) == 0:
                return
            self.edit_peaks = { 'x': edit_x, 'y': edit_y }
        except:
            pass

        try: # mnova
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
                self.edit_peaks = { 'x': edit_x, 'y': edit_y }
        except:
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
        self.edit_peaks = { 'x': edit_x, 'y': edit_y }


    def __run_auto_pick_peak(self):
        max_y = np.max(self.ys)
        height = self.threshold * max_y

        corr_data_ys = self.ys
        corr_height = height
        if 'INFRARED SPECTRUM' == self.datatype:
            corr_data_ys = 1 - self.ys
            corr_height = 1 - height

        peaks = signal.find_peaks(corr_data_ys, height=corr_height)
        peak_idxs = peaks[0]
        self.__set_auto_peaks(peak_idxs)


    def __set_datatable(self):
        y_factor = self.factor and self.factor['y']
        y_factor = y_factor or 1.0
        return encode_datatable(
            self.ys,
            self.boundary['y']['max'],
            y_factor
        )


    def __read_peak_from_file(self):
        self.__read_auto_peaks()
        self.__read_edit_peaks()
        if not self.auto_peaks or not self.params['delta'] == 0.0:
            self.__run_auto_pick_peak()
        if self.params['peaks_str']:
            self.__parse_edit()
