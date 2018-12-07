import nmrglue as ng
import numpy as np
from scipy import signal


THRESHOLD_IR = 0.8
THRESHOLD_NMR = 0.02


class SpectraPeaker():
    def __init__(self, path):
        self.dic, self.data = self.__read(path)
        self.datatype = self.dic['DATATYPE']
        self.title = self.dic['TITLE'][0]
        self.target_idx = self.__index_target()
        self.block_count = self.__count_block()
        self.threshold = self.__thres()
        self.obs_freq = self.__set_obs_freq()
        self.x_unit = self.__set_x_unit()
        self.y = self.__read_y()
        self.x = self.__read_x()
        self.boundary = self.__find_boundary()
        self.label = self.__set_label()
        self.auto_peaks = None
        self.edit_peaks = None


    def __read(self, path):
        return ng.jcampdx.read(path, show_all_data=True, read_err='ignore')


    def __thres(self):
        dts = self.datatype
        if 'NMR SPECTRUM' in dts:
            return THRESHOLD_NMR
        elif 'INFRARED SPECTRUM' in dts:
            return THRESHOLD_IR
        return 0.5


    def __index_target(self):
        target_topics = ['NMR SPECTRUM', 'INFRARED SPECTRUM']
        for tp in target_topics:
            if tp in self.datatype:
                idx = self.datatype.index(tp)

        if 'LINK' in self.datatype:
            idx -= 1

        return idx


    def __count_block(self):
        count = 1
        try:
            count = int(self.dic['BLOCKS'][0])
        except:
            pass

        return count


    def __read_x(self): # TBD
        beg_pt = None
        end_pt = None
        idx = self.target_idx

        if beg_pt is None:
            try:
                beg_pt = float(self.dic['FIRSTX'][idx])
                end_pt = float(self.dic['LASTX'][idx])
            except:
                pass

        if beg_pt is None:
            try:
                obs_freq = float(self.dic['$SFO1'][idx])
                shift = float(self.dic['$OFFSET'][idx])
                beg_pt = float(self.dic['FIRST'][idx].replace(' ', '').split(',')[0]) / obs_freq
                end_pt = float(self.dic['LAST'][idx].replace(' ', '').split(',')[0]) / obs_freq
                shift = beg_pt - shift
                beg_pt = beg_pt - shift
                end_pt = end_pt - shift
            except:
                pass

        num_pt = self.y.shape[0]
        x = np.linspace(beg_pt, end_pt, num=num_pt, endpoint=True)

        if self.x_unit == 'HZ':
            x = x / self.obs_freq

        return x


    def __read_y(self):
        y = None
        if self.target_idx > 0:
            y = self.data['real'][self.target_idx]
        else:
            y = self.data[:]

        return y


    def __find_boundary(self):
        return {
            'x': {
                'max': self.x.max(),
                'min': self.x.min(),
            },
            'y': {
                'max': self.y.max(),
                'min': self.y.min(),
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
            auto_x.append(self.x[idx])
            auto_y.append(self.y[idx])
        if len(auto_x) == 0:
            return
        self.auto_peaks = { 'x': auto_x, 'y': auto_y }


    def __read_auto_peaks(self):
        try:
            auto_x = []
            auto_y = []
            pas_list = self.dic['PEAKASSIGNMENTS']
            pas_idx = len(pas_list) - 2
            pas_idx = 0 if pas_idx < 0 else pas_idx
            pas_target = pas_list[pas_idx]

            pas = pas_target.split('\n')[1:]
            for pa in pas:
                info = pa.replace('(', '').replace(')', '').replace(' ', '').split(',')
                auto_x.append(float(info[0]))
                auto_y.append(float(info[1]))
            if len(auto_x) == 0:
                return
            self.auto_peaks = { 'x': auto_x, 'y': auto_y }
        except:
            pass


    def __read_edit_peaks(self):
        try:
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


    def __parse_edit(self, peaks_str):
        edit_x = []
        edit_y = []
        for p in peaks_str.split('#'):
            info = p.split(',')
            edit_x.append(float(info[0]))
            edit_y.append(float(info[1]))
        if len(edit_x) == 0:
            return
        self.edit_peaks = { 'x': edit_x, 'y': edit_y }


    def __run_pick_peak(self):
        self.__read_auto_peaks()
        self.__read_edit_peaks()
        max_y = np.max(self.y)
        height = self.threshold * max_y

        corr_data_ys = self.y
        corr_height = height
        if 'INFRARED SPECTRUM' in self.datatype:
            corr_data_ys = 1 - self.y
            corr_height = 1 - height

        peaks = signal.find_peaks(corr_data_ys, height=corr_height)
        peak_idxs = peaks[0]
        self.__set_auto_peaks(peak_idxs)


    def pick_peak(self):
        self.__read_auto_peaks()
        self.__read_edit_peaks()

        if (not self.auto_peaks) and (not self.edit_peaks):
            self.__run_pick_peak()


    def write_edit_peak(self, peaks_str):
        err = True
        try:
            self.__read_auto_peaks()
            self.__parse_edit(peaks_str)
            err = False
        except:
            pass

        if (not self.auto_peaks) and (not self.edit_peaks):
            err = True

        return err
