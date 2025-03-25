import hashlib
import subprocess as sbp
import time
import shutil
import numpy as np
import pymzml
import pyopenms

from pathlib import Path
from datetime import datetime

from chem_spectra.lib.shared.buffer import store_byte_in_tmp
from chem_spectra.lib.converter.share import reduce_pts

MARGIN = 1

tmp_dir = Path('./chem_spectra/tmp')  # TBD


class MSConverter:
    def __init__(self, file, params=False):
        self.exact_mz, self.edit_scan, self.thres, param_ext = self.__set_params(params)  # noqa
        self.bound_high = self.exact_mz + MARGIN
        self.bound_low = self.exact_mz - MARGIN
        self.typ = 'MS'
        self.dic = {}
        # - - - - - - - - - - -
        fn = file.name.split('.')
        self.fname = fn[0]
        self.ext = param_ext or fn[-1].lower()
        self.target_dir, self.hash_str = self.__mk_dir()
        self.__get_mzml(file)
        self.runs, self.spectra, self.auto_scan = self.__read_mz_ml()
        self.datatables = self.__set_datatables()
        self.__clean()

    def __set_params(self, params):
        exact_mz = params.get('mass', 0) if params else 0
        edit_scan = params.get('scan', 0) if params else 0
        thres = (params and params.get('thres', 5)) or 5
        ext = params.get('ext', '') if params else ''
        return exact_mz, edit_scan, thres, ext

    def __get_mzml(self, file):
        b_content = file.bcore
        if self.ext == 'raw':
            self.tf = store_byte_in_tmp(
                b_content,
                prefix=self.fname,
                suffix='.RAW',
                directory=self.target_dir.absolute().as_posix()
            )
            self.cmd_msconvert = self.__build_cmd_msconvert()
            self.__run_cmd()
        elif self.ext == 'mzml':
            self.tf = store_byte_in_tmp(
                b_content,
                prefix=self.fname,
                suffix='.mzML',
                directory=self.target_dir.absolute().as_posix()
            )
        elif self.ext == 'mzxml':
            self.tf = store_byte_in_tmp(
                b_content,
                prefix=self.fname,
                suffix='.mzXML',
                directory=self.target_dir.absolute().as_posix()
            )
            exp = pyopenms.MSExperiment()
            pyopenms.MzXMLFile().load(self.tf.name, exp)
            target_path = self.__get_mzml_path().absolute().as_posix()
            pyopenms.MzMLFile().store(target_path, exp)

    def __mk_dir(self):
        hash_str = '{}{}'.format(datetime.now(), self.fname)
        hash_str = str.encode(hash_str)
        hash_str = hashlib.md5(hash_str).hexdigest()
        target_dir = tmp_dir / hash_str
        target_dir.mkdir(parents=True, exist_ok=True)
        return target_dir, hash_str

    def __build_cmd_msconvert(self):
        cmd_msconvert = [
            'docker',
            'exec',
            '-d',
            'msconvert_docker',
            'wine',
            'msconvert',
            '/data/{}/{}'.format(self.hash_str, self.tf.name.split('/')[-1]),
            '-o',
            '/data/{}'.format(self.hash_str),
            '--simAsSpectra',
            '--32',
            '--zlib',
            '--filter',
            '"peakPicking true 1-"',
            '--filter',
            '"zeroSamples removeExtra"',
            '--ignoreUnknownInstrumentError'
        ]
        return cmd_msconvert

    def __run_cmd(self):
        sbp.run(self.cmd_msconvert, timeout=10)

    def __get_mzml_path(self):
        fname = ''.join(self.tf.name.split('/')[-1])
        fname = '.'.join(fname.split('.')[:-1])
        fname = '{}.mzML'.format(fname)
        mzml_path = self.target_dir / fname
        return mzml_path

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

        return ratio, noise_ratio, max_seed
    
    def __get_best_ratio(self, old_ratio, new_ratio, noise_ratio, current_index, current_backup_idx, curr_backup_ratio, old_y, new_y):
        best_ratio, best_idx, backup_ratio, backup_idx = old_ratio, current_index, curr_backup_ratio, current_backup_idx
        best_y = old_y
        if (best_ratio < new_ratio) and (noise_ratio <= 50.0):
            best_idx = current_index
            best_ratio = new_ratio
            best_y = new_y
        elif (new_ratio == 100.0) and (noise_ratio <= 50.0) and (best_y < new_y):
            best_idx = current_index
            best_ratio = new_ratio
            best_y = new_y

        if (backup_ratio < new_ratio):
            backup_idx = current_index
            backup_ratio = new_ratio
        
        return best_ratio, best_idx, backup_ratio, backup_idx, best_y


    def __decode(self, runs, decoded_count=1):
        spectra = []
        best_ratio, best_idx, backup_ratio, backup_idx = 0, 0, 0, 0
        best_y = 0
        # print('this')
        # for idx, data in enumerate(runs):
        #     try:
        #         spc = data.peaks('raw')
        #         spectra.append(reduce_pts(spc))
        #     except:
        #         spectra.append(np.array([]))
        #         continue
            

        #     ratio, noise_ratio, y = self.__get_ratio(spc)
        #     best_ratio, best_idx, backup_ratio, backup_idx, best_y = self.__get_best_ratio(
        #         old_ratio=best_ratio,
        #         new_ratio=ratio,
        #         noise_ratio=noise_ratio,
        #         current_index=idx,
        #         current_backup_idx=backup_idx,
        #         curr_backup_ratio=backup_ratio,
        #         old_y=best_y,
        #         new_y=y
        #     )
        # print('this2')


        if decoded_count == 1:
            for idx, data in enumerate(runs):
                try:
                    spc = data.peaks('raw')
                    spectra.append(reduce_pts(spc))
                except:
                    spectra.append(np.array([]))
                    continue

                

                ratio, noise_ratio, y = self.__get_ratio(spc)
                best_ratio, best_idx, backup_ratio, backup_idx, best_y = self.__get_best_ratio(
                    old_ratio=best_ratio,
                    new_ratio=ratio,
                    noise_ratio=noise_ratio,
                    current_index=idx,
                    current_backup_idx=backup_idx,
                    curr_backup_ratio=backup_ratio,
                    old_y=best_y,
                    new_y=y
                )
        else:
            spectrum_count = runs.get_spectrum_count()
            for idx in range(spectrum_count):
                try:
                    data = runs[idx+1]
                except Exception as e:
                    # cannot retrieve data from scan id, just add an empty spectra
                    spectra.append(np.array([]))
                    continue

                spc = data.peaks('raw')
                spectra.append(reduce_pts(spc))

                ratio, noise_ratio, y = self.__get_ratio(spc)
                # if (best_ratio < ratio) and (noise_ratio <= 50.0):
                #     best_idx = idx
                #     best_ratio = ratio
                #     best_y = y
                # elif (ratio == 100.0) and (noise_ratio <= 50.0) and (best_y < y):
                #     best_idx = idx
                #     best_ratio = ratio
                #     best_y = y

                # if (backup_ratio < ratio):
                #     backup_idx = idx
                #     backup_ratio = ratio
                best_ratio, best_idx, backup_ratio, backup_idx, best_y = self.__get_best_ratio(
                    old_ratio=best_ratio,
                    new_ratio=ratio,
                    noise_ratio=noise_ratio,
                    current_index=idx,
                    current_backup_idx=backup_idx,
                    curr_backup_ratio=backup_ratio,
                    old_y=best_y,
                    new_y=y
                )

        output_idx = best_idx if best_ratio > 10.0 else backup_idx

        return spectra, (output_idx + 1)

    def __read_mz_ml(self):
        mzml_path = self.__get_mzml_path()
        mzml_file = mzml_path.absolute().as_posix()

        runs, spectra, auto_scan = None, None, 0
        elapsed = 0.0
        decoded_count = 1
        while True:
            if mzml_path.exists():
                try:
                    elapsed += 0.2
                    time.sleep(0.2)
                    runs = pymzml.run.Reader(mzml_file, build_index_from_scratch=True)
                    spectra, auto_scan = self.__decode(runs, decoded_count)
                    break
                except:  # noqa
                    decoded_count += 1
                    pass
            else:
                elapsed += 0.1
                time.sleep(0.1)
            if elapsed > 120.0:
                break

        return runs, spectra, auto_scan

    def __set_datatables(self):
        dts = []
        for idx, spc in enumerate(self.spectra):
            # RESOLVE_VSMBNAN2 a valid spectrum must be np.array (N, 2)
            if not spc.shape[0] > 0:
                spc = np.array([[1000.0, 0.0], [2000.0, 0.0]])  # placeholder
                if self.auto_scan == (idx + 1) and (idx < len(self.spectra) - 1):  # move selected scan
                    self.auto_scan += 1
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

    def __clean(self):
        self.tf.close()
        shutil.rmtree(self.target_dir.absolute().as_posix())
