import json
import zipfile
import tempfile
import glob     # noqa: F401
import os
from pathlib import Path

from chem_spectra.lib.shared.buffer import store_str_in_tmp, store_byte_in_tmp
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter
from chem_spectra.lib.converter.jcamp.ms import JcampMSConverter
from chem_spectra.lib.converter.cdf.base import CdfBaseConverter
from chem_spectra.lib.converter.cdf.ms import CdfMSConverter
from chem_spectra.lib.converter.fid.base import FidBaseConverter
from chem_spectra.lib.converter.fid.bruker import FidHasBruckerProcessed
from chem_spectra.lib.converter.bagit.base import BagItBaseConverter
from chem_spectra.lib.converter.lcms.lcms import LCMSConverter
from chem_spectra.lib.converter.ms import MSConverter
from chem_spectra.lib.composer.ni import NIComposer
from chem_spectra.lib.composer.ms import MSComposer
from chem_spectra.lib.composer.lcms import LCMSComposer
from chem_spectra.lib.composer.base import BaseComposer     # noqa: F401
from chem_spectra.lib.converter.nmrium.base import NMRiumDataConverter
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.path as mpath  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from chem_spectra.model.concern.property import decorate_sim_property
from chem_spectra.lib.external.binaryparser import get_openlab_readers
from chem_spectra.lib.external.chemotion_converter_lcms import lcms_frames_from_converter_app


def find_dir(path, name):
    for root, _, files in os.walk(path):
        if name in files:
            return os.path.join(root)
    return False


def search_brucker_binary(td):
    try:
        has_processed_files = search_processed_file(td)
        target_dir = find_dir(td, 'fid')
        if target_dir:
            return target_dir, has_processed_files
        return False, has_processed_files
    except Exception:  # noqa: E722
        return False, False

def search_processed_file(td):
    try:
        pdata_dir = find_and_get_dir(td, 'pdata')
        for root, dirs, _ in os.walk(pdata_dir):
            return (len(dirs) > 0)
    except:
        return False

def find_and_get_dir(path, name):
    for root, dirs, _ in os.walk(path):
        if name in dirs:
            return os.path.join(root, name)
    return False

def search_bag_it_file(td):
    try:
        target_dir = find_dir(td, 'bagit.txt')
        return target_dir
    except:     # noqa: E722
        return False

def search_lcms_file(td):
    try:
        target_dir = find_dir(td, 'LCMS.csv')
        return target_dir
    except:     # noqa: E722
        return False

def _find_dir_with_cdf(root: str):
    for dirpath, _, filenames in os.walk(root):
        for fn in filenames:
            if fn.lower().endswith('.cdf'):
                return dirpath
    return None

def compute_tic(ms_df: pd.DataFrame) -> pd.DataFrame:
    return (
        ms_df.groupby('time', as_index=False)['intensities'].sum()
        .sort_values('time')
        .reset_index(drop=True)
        .rename(columns={'intensities': 'Intensity'})
    )

class TransformerModel:
    def __init__(self, file, molfile=None, params=False, multiple_files=False):
        self.file = file
        self.molfile = molfile
        self.params = params
        self.multiple_files = multiple_files

    def convert2jcamp(self):
        cmpsr, _ = self.to_composer()
        if isinstance(cmpsr, BagItBaseConverter):
            return cmpsr
        return cmpsr.tf_jcamp()

    def convert2img(self):
        cmpsr, _ = self.to_composer()
        if isinstance(cmpsr, BagItBaseConverter):
            return cmpsr
        return cmpsr.tf_img()

    def convert2jcamp_img(self):
        cmpsr, _ = self.to_composer()
        if not cmpsr:
            return False, False, False
        if isinstance(cmpsr, BagItBaseConverter):
            return False, cmpsr, False
        return cmpsr.tf_jcamp(), cmpsr.tf_img(), cmpsr.tf_csv()

    def to_composer(self):
        is_raw_mzml = self.file.name.split('.')[-1].lower() in ['raw', 'mzml', 'mzxml']     # noqa: E501
        is_cdf = self.file.name.split('.')[-1].lower() in ['cdf']
        is_zip = self.file.name.split('.')[-1].lower() in ['zip']
        is_raw_mzml_by_params = self.params['ext'] in ['raw', 'mzml', 'mzxml']
        is_cdf_by_params = self.params['ext'] in ['cdf']
        is_zip_by_params = self.params['ext'] in ['zip']
        if is_raw_mzml or is_raw_mzml_by_params:
            return self.ms2composer(), False
        if is_cdf or is_cdf_by_params:
            _, cp = self.cdf2cvp()
            return cp, False
        if is_zip or is_zip_by_params:
            _, cp, invalid_molfile = self.zip2cvp()
            return cp, invalid_molfile
        else:
            _, cp, invalid_molfile = self.jcamp2cvp()
            return cp, invalid_molfile

    def to_converter(self):
        is_raw_mzml = self.file.name.split('.')[-1].lower() in ['raw', 'mzml', 'mzxml']     # noqa: E501
        is_cdf = self.file.name.split('.')[-1].lower() in ['cdf']
        is_zip = self.file.name.split('.')[-1].lower() in ['zip']
        is_raw_mzml_by_params = self.params['ext'] in ['raw', 'mzml', 'mzxml']
        is_cdf_by_params = self.params['ext'] in ['cdf']
        is_zip_by_params = self.params['ext'] in ['zip']
        if is_raw_mzml or is_raw_mzml_by_params:
            return self.ms2composer()
        if is_cdf or is_cdf_by_params:
            cv, _ = self.cdf2cvp()
            return cv
        if is_zip or is_zip_by_params:
            cv, _, _ = self.zip2cvp()
            return cv
        else:
            cv, _ = self.jcamp2cvp()
            return cv

    def ms2composer(self):
        mscv = MSConverter(self.file, self.params)
        mscp = MSComposer(mscv)
        return mscp

    def cdf2cvp(self):
        tf = store_byte_in_tmp(self.file.bcore, suffix='.cdf')
        cbcv = CdfBaseConverter(tf.name, self.params)
        tf.close()
        # conversion
        mscv = CdfMSConverter(cbcv)
        mscp = MSComposer(mscv)
        return mscv, mscp

    def zip2cvp(self):
        with tempfile.TemporaryDirectory() as td:
            tz = store_byte_in_tmp(self.file.bcore, suffix='.zip')
            with zipfile.ZipFile(tz.name, 'r') as z:
                z.extractall(td)
            target_dir, has_processed_files = search_brucker_binary(td)
            if target_dir:
                # NMR data
                if (has_processed_files):
                    return self.zip2cv_with_processed_file(target_dir, self.params, self.file.name)
                fbcv = FidBaseConverter(target_dir, self.params, self.file.name)
                if not fbcv:
                    return False, False, False

                isSimulateNMR = False
                if self.params and 'simulatenmr' in self.params:
                    isSimulateNMR = self.params['simulatenmr']
                decorated_jbcv = decorate_sim_property(fbcv, self.molfile, isSimulateNMR)   # noqa: E501
                invalid_molfile = False
                if self.molfile is None:
                    invalid_molfile = True
                if ((type(decorated_jbcv) is dict) and "invalid_molfile" in decorated_jbcv):
                    invalid_molfile = True
                    final_decorated_jbcv = decorated_jbcv['origin_jbcv']
                else:
                    final_decorated_jbcv = decorated_jbcv
                nicv = JcampNIConverter(final_decorated_jbcv)
                nicp = NIComposer(nicv)
                return nicv, nicp, invalid_molfile

            openlab_dir = _find_dir_with_cdf(td)
            if openlab_dir:
                normalized_dir = os.path.join(td, 'normalized')
                os.makedirs(normalized_dir, exist_ok=True)

                converter_frames = lcms_frames_from_converter_app(openlab_dir)
                if converter_frames is not None:
                    lc_df, minus_df, plus_df = converter_frames
                else:
                    read_lc, read_ms = get_openlab_readers()
                    if read_lc is None or read_ms is None:
                        raise RuntimeError('OpenLab readers are not available')
                    lc_df = read_lc(openlab_dir)
                    minus_df, plus_df = read_ms(openlab_dir)

                required_lc_cols = ['RetentionTime', 'DetectorSignal', 'wavelength']
                for col in required_lc_cols:
                    if col not in lc_df.columns:
                        raise RuntimeError(f'LC output missing required column {col}')
                lc_df = lc_df[required_lc_cols]

                for label, df in {'MINUS': minus_df, 'PLUS': plus_df}.items():
                    for col in ('mz', 'intensities', 'time'):
                        if col not in df.columns:
                            raise RuntimeError(f'MS {label} missing column {col}')

                tic_minus = compute_tic(minus_df)
                tic_plus = compute_tic(plus_df)

                lc_df.to_csv(os.path.join(normalized_dir, 'LCMS.csv'), index=False)
                tic_minus.to_csv(os.path.join(normalized_dir, 'TIC_MINUS.csv'), index=False)
                tic_plus.to_csv(os.path.join(normalized_dir, 'TIC_PLUS.csv'), index=False)
                minus_df[['time', 'mz', 'intensities']].to_csv(
                    os.path.join(normalized_dir, 'MZ_MINUS_Spectra.csv'), index=False
                )
                plus_df[['time', 'mz', 'intensities']].to_csv(
                    os.path.join(normalized_dir, 'MZ_PLUS_Spectra.csv'), index=False
                )

                lcms_cv = LCMSConverter(normalized_dir, self.params, os.path.basename(self.file.name))
                lcms_peaks = None
                if self.params and 'lcms_peaks' in self.params:
                    try:
                        lcms_peaks = json.loads(self.params['lcms_peaks'])
                    except (json.JSONDecodeError, TypeError):
                        lcms_peaks = None
                lcms_np = LCMSComposer(lcms_cv, lcms_peaks)
                return lcms_cv, lcms_np, False

            is_bagit = search_bag_it_file(td)
            if is_bagit:
                bagcv = BagItBaseConverter(td, self.params, self.file.name)
                return bagcv, bagcv, False

            is_lcms = search_lcms_file(td)
            if is_lcms:
                lcms_path = is_lcms
                lcms_cv = LCMSConverter(lcms_path, self.params, os.path.basename(lcms_path))

                lcms_peaks = None
                if self.params and 'lcms_peaks' in self.params:
                    try:
                        lcms_peaks = json.loads(self.params['lcms_peaks'])
                    except (json.JSONDecodeError, TypeError):
                        lcms_peaks = None

                lcms_np = LCMSComposer(lcms_cv, lcms_peaks)
                return lcms_cv, lcms_np, False

        return False, False, False

    def zip2cv_with_processed_file(self, target_dir, params, file_name):
        fid_brucker = FidHasBruckerProcessed(target_dir, params, file_name)
        if not fid_brucker:
            return False, False, False

        isSimulateNMR = False
        if params and 'simulatenmr' in params:
            isSimulateNMR = params['simulatenmr']

            
        list_decorated_converters = []
        list_decorated_composers = []

        invalid_molfile = False
        for conv in fid_brucker.data:
            decorated_jbcv = decorate_sim_property(conv, self.molfile, isSimulateNMR)   # noqa: E501
            
            if ((type(decorated_jbcv) is dict) and "invalid_molfile" in decorated_jbcv):
                invalid_molfile = True
                final_decorated_jbcv = decorated_jbcv['origin_jbcv']
            else:
                final_decorated_jbcv = decorated_jbcv
            

        # invalid_molfile = False
        # for conv in fid_brucker.data:
        #     if ((type(decorated_jbcv) is dict) and "invalid_molfile" in decorated_jbcv):
        #         invalid_molfile = True
        #         final_decorated_jbcv = decorated_jbcv['origin_jbcv']
        #     else:
        #         final_decorated_jbcv = decorated_jbcv
        #         final_decorated_jbcv.simu_peaks = decorated_jbcv.simu_peaks

            list_decorated_converters.append(final_decorated_jbcv)
            nicv = JcampNIConverter(final_decorated_jbcv)
            nicp = NIComposer(nicv)
            list_decorated_composers.append(nicp)
        return list_decorated_converters, list_decorated_composers, invalid_molfile

    def jcamp2cvp(self):
        tf = store_str_in_tmp(self.file.core)
        jbcv = JcampBaseConverter(tf.name, self.params)
        tf.close()
        invalid_molfile = False
        # conversion
        if jbcv.typ == 'MS':
            mscv = JcampMSConverter(jbcv)
            mscp = MSComposer(mscv)
            return mscv, mscp, invalid_molfile
        elif jbcv.typ == 'LC/MS':
            lcmscv = LCMSConverter(jbcv, params=jbcv.params, fname=(jbcv.fname or 'chromatogram.jdx'))
            lcmscp = LCMSComposer(lcmscv)
            return lcmscv, lcmscp, invalid_molfile
        else:
            isSimulateNMR = False
            if self.params and 'simulatenmr' in self.params:
                isSimulateNMR = self.params['simulatenmr']
            decorated_jbcv = decorate_sim_property(jbcv, self.molfile, isSimulateNMR)   # noqa: E501
            
            if ((type(decorated_jbcv) is dict) and "invalid_molfile" in decorated_jbcv):
                invalid_molfile = True
                final_decorated_jbcv = decorated_jbcv['origin_jbcv']
            else:
                final_decorated_jbcv = decorated_jbcv

            nicv = JcampNIConverter(final_decorated_jbcv)
            nicp = NIComposer(nicv)
            return nicv, nicp, invalid_molfile

    def tf_predict(self):
        target = json.loads(self.params['predict'])
        r = target.get('output', {}).get('result')
        if not r or not r[0]:
            return False

        tf = store_str_in_tmp(
            json.dumps(target),
            suffix='.json',
        )
        return tf

    def tf_nmrium(self):
        converter = NMRiumDataConverter(self.file)
        if converter.is_2d == True or converter.data is None:
            return None
        nicp = NIComposer(converter)
        tf_jcamp = nicp.tf_jcamp()
        return tf_jcamp

    def __get_cyclic_volta_ref_peaks(self, curve_idx, extraParams):
        x_peaks, y_peaks = [], []
        try:
            extras_dict = json.loads(extraParams) if extraParams else None
            cyclic_volta_str = extras_dict['cyclicvolta'] if extras_dict else None
            cyclic_volta = json.loads(cyclic_volta_str) if cyclic_volta_str else None
            spectra_list = cyclic_volta['spectraList'] if cyclic_volta else None
            spectra_extra = spectra_list[curve_idx] if spectra_list and curve_idx < len(spectra_list) else None
            list_peaks = spectra_extra['list'] if spectra_extra else []
            x_peaks, y_peaks = [], []
            for peak in list_peaks:
                min_peak, max_peak, isRef = peak['min'], peak['max'], peak['isRef']
                if isRef == True:
                    x_peaks.extend([min_peak['x'], max_peak['x']])
                    y_peaks.extend([min_peak['y'], max_peak['y']])
        except:
            pass

        return x_peaks, y_peaks

    def tf_combine(self, list_file_names=None, extraParams=None):
        if not self.multiple_files:
            return False

        path_data = [
            (mpath.Path.MOVETO, (0, 5)),
            (mpath.Path.LINETO, (0, 20)),
        ]
        codes, verts = zip(*path_data)

        circle = mpath.Path.unit_circle()
        cirle_verts = np.concatenate([circle.vertices, verts])
        cirle_codes = np.concatenate([circle.codes, codes])
        cut_star_marker = mpath.Path(cirle_verts, cirle_codes)
          
        plt.rcParams['figure.figsize'] = [16, 9]
        plt.rcParams['font.size'] = 14
        plt.rcParams['legend.loc'] = 'upper left'
        curve_idx = self.params.get('jcamp_idx', 0)

        xlabel, ylabel = '', ''
        xlabel_set, ylabel_set = [], []
        dic_x_label, dic_y_label = {}, {}

        for idx, file in enumerate(self.multiple_files):
            if (list_file_names is not None) and idx < len(list_file_names):
                file.name = list_file_names[idx]
                self.multiple_files[idx] = file

        self.multiple_files.sort(key=lambda file: file.name)
        
        for idx, file in enumerate(self.multiple_files):
            tf = store_str_in_tmp(file.core)
            jbcv = JcampBaseConverter(tf.name, self.params)
            filename = file.name
            if jbcv.typ == 'MS':
                mscv = JcampMSConverter(jbcv)
                mscp = MSComposer(mscv)
                plt.plot(mscp.core.xs, mscp.core.ys, label=filename)
            else:
                nicv = JcampNIConverter(jbcv)
                nicp = NIComposer(nicv)
                xs, ys = nicp.core.xs, nicp.core.ys
                marker = ''
                if nicp.core.is_aif:
                    first_x, last_x = xs[0], xs[len(xs)-1]
                    if first_x <= last_x:
                        filename = 'ADSORPTION'
                        marker = '^'
                    else:
                        filename = 'DESORPTION'
                        marker = 'v'
                plt.plot(xs, ys, label=filename, marker=marker)

                # PLOT label
                core_label_x = nicp.core.label['x']
                core_label_y = nicp.core.label['y']
                if nicp.core.is_cyclic_volta:
                    x_peaks, y_peaks = self.__get_cyclic_volta_ref_peaks(curve_idx, extraParams)

                    plt.plot(
                        x_peaks,
                        y_peaks,
                        'r',
                        ls='',
                        marker=cut_star_marker,
                        markersize=50,
                    )

                    if core_label_x not in dic_x_label:
                        xlabel_set.append(core_label_x)
                        dic_x_label[core_label_x] = 1
                    if core_label_y not in dic_y_label:
                        ylabel_set.append(core_label_y)
                        dic_y_label[core_label_y] = 1
                    if (idx == len(self.multiple_files) - 1):
                        xlabel = ', '.join(xlabel_set)
                        ylabel = ', '.join(ylabel_set)
                elif (nicp.core.non_nmr == False):
                    xlabel = "Chemical shift ({})".format(core_label_x.lower())
                    ylabel = "Intensity ({})".format(core_label_y.lower())
                else:
                    xlabel = "X ({})".format(core_label_x)
                    ylabel = "Y ({})".format(core_label_y)

            tf.close()
        
        plt.xlabel(xlabel, fontsize=18)
        plt.ylabel(ylabel, fontsize=18)
        plt.legend()
        tf_img = tempfile.NamedTemporaryFile(suffix='.png')
        plt.savefig(tf_img, format='png')
        tf_img.seek(0)
        plt.clf()
        plt.cla()
        return tf_img
