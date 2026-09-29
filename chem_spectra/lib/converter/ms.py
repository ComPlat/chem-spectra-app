from chem_spectra.lib.converter.jcamp.techniques import technique_for
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

# How long to wait for the msconvert sidecar to answer. It must exceed the
# ceiling the sidecar imposes on itself -- `mscrunner.py` runs the conversion
# under `timeout=10` -- so that its own reply arrives instead of being cut off
# here. When this was also 10s the two raced, and ours usually won, which threw
# away the sidecar's account of what went wrong.
SIDECAR_TIMEOUT = 30

# How long to wait for the .mzML to appear after the sidecar has answered.
#
# In the deployed topology this should be ~0: the shim blocks on its HTTP
# request and `mscrunner` runs msconvert synchronously, so the file exists by
# the time __run_cmd returns. It stays generous because INSTALL.md documents a
# different local setup -- a real `docker exec -d` against the plain
# ProteoWizard image -- where `-d` genuinely detaches and the poll is the only
# thing waiting for the conversion.
MZML_WAIT = 120.0


class MSConversionFailed(RuntimeError):
    """msconvert produced no usable mzML for this upload.

    Raised instead of letting the failure surface further down as
    `TypeError: 'NoneType' object is not iterable`, which said nothing about
    the cause and only appeared after the full MZML_WAIT.

    `status` is what the request answers with. 422 is the refusal convention:
    the upload was understood and cannot be processed. A converter that is not
    answering is a different thing and says so below -- telling a chemist their
    file is unprocessable when the service is down sends them to look in the
    wrong place.
    """

    status = 422


class MSConverterUnavailable(MSConversionFailed):
    """The msconvert sidecar could not be reached, or did not answer in time.

    Nothing is wrong with the upload, so this is ours rather than the
    caller's. It is still delivered as JSON, because that is the only way the
    ELN shows a reason at all.
    """

    status = 502


class MSConverter:
    def __init__(self, file, params=False):
        self.exact_mz, self.edit_scan, self.thres, param_ext = self.__set_params(params)  # noqa
        self.bound_high = self.exact_mz + MARGIN
        self.bound_low = self.exact_mz - MARGIN
        self.typ = 'MS'
        # RAW / mzML / mzXML never go through JCAMP classification; typ is
        # stated literally above, so the descriptor is known here.
        self.technique = technique_for(self.typ)
        self.dic = {}
        # - - - - - - - - - - -
        fn = file.name.split('.')
        self.fname = fn[0]
        self.ext = param_ext or fn[-1].lower()
        self.target_dir, self.hash_str = self.__mk_dir()
        # __clean in a finally: it used to be the last statement here, so any
        # failure above left the uploaded file and its hashed directory on the
        # shared volume for good.
        try:
            self.__get_mzml(file)
            self.runs, self.spectra, self.auto_scan = self.__read_mz_ml()
            self.datatables = self.__set_datatables()
        finally:
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
        """Hand the command to the sidecar, and keep what it says.

        The result used to be discarded entirely -- return code, stdout and
        stderr. It is worth keeping, because `/bin/docker` in the deployed
        image is not Docker: it is a shim that POSTs this command to the
        msconvert service, and its `requests.post` sits outside its own
        try/except, so it exits non-zero when the service cannot be reached.
        That is precisely the failure that used to appear two minutes later as
        an unexplained `TypeError`.

        The `-d` in the command is inert -- there is no daemon, and the shim
        blocks on its HTTP request -- but it must stay: with MSC_VALIDATE set,
        which the deployed images do, the shim asserts on it.
        """
        try:
            result = sbp.run(
                self.cmd_msconvert,
                timeout=SIDECAR_TIMEOUT,
                capture_output=True,
                text=True,
            )
        except sbp.TimeoutExpired:
            raise MSConverterUnavailable(
                'the msconvert service did not answer within '
                '{}s'.format(SIDECAR_TIMEOUT))

        if result.returncode != 0:
            detail = (result.stderr or result.stdout or '').strip()
            raise MSConverterUnavailable(
                'the msconvert service could not be reached, or refused the '
                'command{}'.format(': ' + detail if detail else ''))

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

        # Only the RAW path has a conversion to wait for. mzML and mzXML are
        # written synchronously by __get_mzml, so a file that will not parse
        # will not parse in two minutes either; waiting only delays the answer
        # and makes it look like a converter problem.
        wait_for = MZML_WAIT if self.ext == 'raw' else 0.0

        runs, spectra, auto_scan = None, None, 0
        elapsed = 0.0
        decoded_count = 1
        unparsable = None
        while True:
            if mzml_path.exists():
                try:
                    elapsed += 0.2
                    time.sleep(0.2)
                    runs = pymzml.run.Reader(mzml_file, build_index_from_scratch=True)
                    spectra, auto_scan = self.__decode(runs, decoded_count)
                    break
                except Exception as err:  # noqa
                    decoded_count += 1
                    unparsable = err
            else:
                elapsed += 0.1
                time.sleep(0.1)
            if elapsed > wait_for:
                raise self.__no_spectra(mzml_path, unparsable)

        return runs, spectra, auto_scan

    def __no_spectra(self, mzml_path, unparsable):
        """Say which of the two failures this is.

        They are different problems for whoever has to act on them: a file
        that never arrived is the converter's, a file that arrived and will
        not parse is the data's. Reporting both as "no mzML appeared" sent the
        reader after the wrong one -- and for an mzML upload, where nothing is
        converted at all, it blamed a converter that never ran.
        """
        if unparsable is not None:
            return MSConversionFailed(
                '{} could not be read as mzML: {}: {}'.format(
                    mzml_path.name, type(unparsable).__name__, unparsable))
        return MSConversionFailed(
            'no mzML was produced for this upload: {} never appeared, '
            '{:.0f}s after msconvert reported success'.format(
                mzml_path.name, MZML_WAIT))

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
        """Tolerant on purpose: this now runs in a finally.

        `self.tf` is only set for the extensions __get_mzml knows, and the
        directory may be half-built, so a strict cleanup would replace the real
        exception with its own.
        """
        handle = getattr(self, 'tf', None)
        if handle is not None:
            handle.close()
        shutil.rmtree(self.target_dir.absolute().as_posix(), ignore_errors=True)
