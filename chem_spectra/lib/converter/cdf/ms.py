import numpy as np

from chem_spectra.lib.converter.jcamp.techniques import technique_for
from chem_spectra.lib.shared.ms_scan import scan_ratios
from chem_spectra.lib.shared.number_text import xy_lines

MARGIN = 1
THRESHOLD_MS = 0.05


class CdfMSConverter:  # nmr & IR
    def __init__(self, base):
        self.params = base.params
        self.dic = base.dic
        self.data = base.data
        self.datatypes = base.datatypes
        self.datatype = base.datatype
        self.title = base.title
        self.typ = base.typ
        # CdfBaseConverter hardcodes typ='MS', so the descriptor is known
        # here without consulting data_type.json -- this path never goes
        # through JCAMP classification at all.
        self.technique = technique_for(self.typ)
        # - - - - - - - - - - -
        self.exact_mz, self.edit_scan, self.thres = self.__set_params(base.params)  # noqa
        self.bound_high = self.exact_mz + MARGIN
        self.bound_low = self.exact_mz - MARGIN
        # - - - - - - - - - - -
        self.fname = base.fname
        self.runs, self.spectra, self.auto_scan = self.__read_data()
        self.datatables = self.__set_datatables()

    def __set_params(self, params):
        exact_mz = params.get('mass', 0) if params else 0
        pscan = params.get('scan', None) if params else None
        escan = self.dic.get('$CSSCANEDITTARGET', [])
        escan = int(escan[0]) if len(escan) > 0 else None
        edit_scan = pscan or escan or None
        pthres = params.get('thres', None) if params else None
        ethres = self.dic.get('$CSTHRESHOLD', [])
        ethres = float(ethres[0]) * 100.0 if len(ethres) > 0 else None
        thres = pthres or ethres or (THRESHOLD_MS * 100.0)
        return exact_mz, edit_scan, thres

    def __get_ratio(self, spc):
        ratio, noise_ratio, max_seed = scan_ratios(
            spc, self.bound_low, self.bound_high)
        return ratio, noise_ratio

    def __decode(self, runs):
        spectra = []
        best_ratio, best_idx, backup_ratio, backup_idx = 0, 0, 0, 0
        for idx, data in enumerate(runs):
            spectra.append(data)

            ratio, noise_ratio = self.__get_ratio(data)
            if (best_ratio < ratio) and (noise_ratio <= 50.0):
                best_idx = idx
                best_ratio = ratio

            if (backup_ratio < ratio):
                backup_idx = idx
                backup_ratio = ratio

        output_idx = best_idx if best_ratio > 10.0 else backup_idx

        return spectra, (output_idx + 1)

    def __build_runs(self):
        data = self.data
        runs = []
        xs, ys, idxs = data['XS'], data['YS'], data['SCANINDEXS']
        for i, v in enumerate(idxs[:-1]):
            head, tail = v, idxs[i+1]
            dxs, dys = xs[head:tail].data, ys[head:tail].data
            ds = np.array([dxs, dys])
            ds = np.transpose(ds)
            runs.append(ds)
        return runs

    def __read_data(self):
        runs, spectra, auto_scan = self.__build_runs(), None, 0
        spectra, auto_scan = self.__decode(runs)
        return runs, spectra, auto_scan

    def __set_datatables(self):
        dts = []
        for idx, spc in enumerate(self.spectra):
            xs = spc[:, 0]
            ys = spc[:, 1]
            pts = xs.shape[0]
            dt = xy_lines(xs, ys)
            dts.append({'dt': dt, 'pts': pts})
        return dts
