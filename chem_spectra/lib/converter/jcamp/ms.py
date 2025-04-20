from chem_spectra.lib.converter.jcamp.data_parse import make_ms_data_xsys
from chem_spectra.lib.converter.share import reduce_pts
import numpy as np

MARGIN = 1
THRESHOLD_MS = 0.05


class JcampMSConverter:  # nmr & IR
    def __init__(self, base):
        self.params = base.params
        self.datatypes = base.datatypes
        self.datatype = base.datatype
        self.typ = base.typ
        self.dic = base.dic
        self.data = make_ms_data_xsys(base)
        self.title = base.title
        self.is_em_wave = base.is_em_wave
        self.is_ir = base.is_ir
        self.non_nmr = base.non_nmr
        # - - - - - - - - - - -
        self.exact_mz, self.edit_scan, self.thres = self.__set_params(base.params)  # noqa
        self.bound_high = self.exact_mz + MARGIN
        self.bound_low = self.exact_mz - MARGIN
        # - - - - - - - - - - -
        self.fname = base.fname
        self.runs, self.spectra, self.auto_scan = self.__read_mz_ml()
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
        all_ys, ratio = [], 0
        bLow, bHigh = self.bound_low, self.bound_high

        match_base_xs, match_base_ys = [], []
        match_seed_xs, match_seed_ys = [], []
        match_oorg_xs, match_oorg_ys = [], []

        for pk in spc:
            all_ys.append(pk[1])

            if bLow < pk[0] < bHigh:
                match_seed_xs.append(pk[0])
                match_seed_ys.append(pk[1])
            elif bLow + 1 < pk[0] < bHigh + 1:
                match_seed_xs.append(pk[0])
                match_seed_ys.append(pk[1])
            elif bLow + 23 < pk[0] < bHigh + 23:
                match_seed_xs.append(pk[0])
                match_seed_ys.append(pk[1])
            elif bLow + 39 < pk[0] < bHigh + 39:
                match_seed_xs.append(pk[0])
                match_seed_ys.append(pk[1])

            if pk[0] <= bHigh + 39:
                match_base_xs.append(pk[0])
                match_base_ys.append(pk[1])
            elif bHigh + 39 < pk[0]:
                match_oorg_xs.append(pk[0])
                match_oorg_ys.append(pk[1])

        max_base = max(match_base_ys, default=0.1)
        max_seed = max(match_seed_ys, default=0)
        max_oorg = max(match_oorg_ys, default=0)

        ratio = 100 * max_seed / max_base
        noise_ratio = 100 * max_oorg / max_base

        return ratio, noise_ratio

    # Guarantees that the data passed to `reduce_pts` and subsequent operations
    # conforms to the expected structure.
    def __normalize_data(self, data):
        if isinstance(data, np.ndarray) and data.ndim == 3 and data.shape[0] == 1:
            return data[0]
        return data

    def __decode(self, runs):
        spectra = []
        best_ratio, best_idx, backup_ratio, backup_idx = 0, 0, 0, 0
        for idx, data in enumerate(runs):
            data = self.__normalize_data(data)
            spectra.append(reduce_pts(data))
            ratio, noise_ratio = self.__get_ratio(data)
            if (best_ratio < ratio) and (noise_ratio <= 50.0):
                best_idx = idx
                best_ratio = ratio

            if (backup_ratio < ratio):
                backup_idx = idx
                backup_ratio = ratio

        output_idx = best_idx if best_ratio > 10.0 else backup_idx

        return spectra, (output_idx + 1)

    def __read_mz_ml(self):
        runs, spectra, auto_scan = self.data, None, 0
        spectra, auto_scan = self.__decode(runs)
        return runs, spectra, auto_scan

    def __set_datatables(self):
        dts = []
        for idx, spc in enumerate(self.spectra):
            xs = spc[:, 0]
            ys = spc[:, 1]
            pts = xs.shape[0]
            dt = []
            for idx in range(pts):
                dt.append(
                    '{}, {}\n'.format(
                        xs[idx],
                        ys[idx]
                    )
                )
            dts.append({'dt': dt, 'pts': pts})
        return dts
