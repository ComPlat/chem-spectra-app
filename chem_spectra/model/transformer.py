import json
import zipfile
import tempfile
import glob     # noqa: F401
import os

from chem_spectra.lib.shared.buffer import store_str_in_tmp, store_byte_in_tmp
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter
from chem_spectra.lib.converter.jcamp.ms import JcampMSConverter
from chem_spectra.lib.converter.cdf.base import CdfBaseConverter
from chem_spectra.lib.converter.cdf.ms import CdfMSConverter
from chem_spectra.lib.converter.fid.base import FidBaseConverter
from chem_spectra.lib.converter.fid.bruker import FidHasBruckerProcessed
from chem_spectra.lib.converter.bagit.base import BagItBaseConverter
from chem_spectra.lib.converter.ms import MSConverter
from chem_spectra.lib.composer.ni import NIComposer
from chem_spectra.lib.composer.ms import MSComposer
from chem_spectra.lib.composer.base import BaseComposer     # noqa: F401
from chem_spectra.lib.converter.nmrium.base import NMRiumDataConverter
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import ticker  # noqa: E402

from chem_spectra.model.concern.property import decorate_sim_property


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
            cv, _ = self.zip2cvp()
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
        fbcv = False
        with tempfile.TemporaryDirectory() as td:
            tz = store_byte_in_tmp(self.file.bcore, suffix='.zip')
            with zipfile.ZipFile(tz.name, 'r') as z:
                z.extractall(td)
            target_dir, has_processed_files = search_brucker_binary(td)
            if target_dir:
                # NMR data
                if (has_processed_files):
                    return self.zip2cv_with_processed_file(target_dir, self.params, self.file.name)
                else:
                    fbcv = FidBaseConverter(target_dir, self.params, self.file.name)
                    if not fbcv:
                        return False, False, False

                    isSimulateNMR = False
                    if self.params and 'simulatenmr' in self.params:
                        isSimulateNMR = self.params['simulatenmr']
                    decorated_jbcv = decorate_sim_property(fbcv, self.molfile, isSimulateNMR)   # noqa: E501

                    # if ((type(decorated_jbcv) is dict) and "invalid_molfile" in decorated_jbcv):
                    #     # return if molfile is invalid
                    #     return None, decorated_jbcv
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
            else:
                is_bagit = search_bag_it_file(td)
                if is_bagit:
                    bagcv = BagItBaseConverter(td, self.params, self.file.name)
                    return bagcv, bagcv, False

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

    def tf_combine(self, list_file_names=None, extraParams=None):
        if not self.multiple_files:
            return False

        plt.rcParams['figure.figsize'] = [16, 9]
        plt.rcParams['figure.dpi'] = 200
        plt.rcParams['font.size'] = 14
        plt.rcParams['legend.loc'] = 'upper left'
        curve_idx = self.params.get('jcamp_idx', 0)

        xlabel, ylabel = '', ''
        xlabel_set, ylabel_set = [], []
        dic_x_label, dic_y_label = {}, {}
        global_x_min, global_x_max = None, None
        any_forward_orientation = False

        cv_mode = False
        cv_abs_max = 0.0
        for idx, file in enumerate(self.multiple_files):
            if (list_file_names is not None) and idx < len(list_file_names):
                file.name = list_file_names[idx]
                self.multiple_files[idx] = file

        active_name = None
        if 0 <= curve_idx < len(self.multiple_files):
            active_name = self.multiple_files[curve_idx].name

        self.multiple_files.sort(key=lambda file: file.name)

        active_nicp = None
        active_y_values = None

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
                y_values = ys
                is_active = active_name is not None and file.name == active_name
                if nicp.core.is_cyclic_volta:
                    cv_state = {}
                    if extraParams:
                        try:
                            extras_dict = json.loads(extraParams)
                        except Exception:
                            extras_dict = {}
                        cv_state = (
                            extras_dict.get('cyclicvoltaSt')
                            or extras_dict.get('cyclicvolta')
                            or extras_dict.get('cyclic_volta')
                        ) or {}
                    if not cv_state:
                        cv_state = (
                            nicp.core.params.get('cyclicvoltaSt')
                            or nicp.core.params.get('cyclicvolta')
                            or nicp.core.params.get('cyclic_volta')
                        ) or {}
                    if isinstance(cv_state, str):
                        try:
                            cv_state = json.loads(cv_state)
                        except Exception:
                            cv_state = {}
                    cv_display = cv_state.get('cvDisplay') or {}
                    if isinstance(cv_display, str):
                        try:
                            cv_display = json.loads(cv_display)
                        except Exception:
                            cv_display = {}
                    try:
                        scale = float(cv_display.get('yScaleFactor', 1.0))
                    except Exception:
                        scale = 1.0
                    if scale != 1.0:
                        y_values = ys * scale
                    cv_mode = True
                    try:
                        cv_abs_max = max(cv_abs_max, float(np.max(np.abs(y_values))))
                    except Exception:
                        pass
                marker = ''
                if nicp.core.is_aif:
                    first_x, last_x = xs[0], xs[len(xs)-1]
                    if first_x <= last_x:
                        filename = 'ADSORPTION'
                        marker = '^'
                    else:
                        filename = 'DESORPTION'
                        marker = 'v'
                plt.plot(xs, y_values, label=filename, marker=marker)

                if is_active:
                    active_nicp = nicp
                    active_y_values = y_values
                    if nicp.core.is_cyclic_volta:
                        nicp._cv_density_scale = scale

                try:
                    x_max = np.max(xs)
                    x_min = np.min(xs)
                    global_x_min = x_min if global_x_min is None else min(global_x_min, x_min)
                    global_x_max = x_max if global_x_max is None else max(global_x_max, x_max)
                except Exception:
                    pass
                if (nicp.core.is_tga or nicp.core.is_gc or nicp.core.is_uv_vis or nicp.core.is_hplc_uv_vis
                    or nicp.core.is_xrd or nicp.core.is_cyclic_volta or nicp.core.is_sec
                    or nicp.core.is_cds or nicp.core.is_aif or nicp.core.is_emissions
                    or nicp.core.is_dls_acf or nicp.core.is_dls_intensity):
                    any_forward_orientation = True

                # PLOT label
                core_label_x = nicp.core.label['x']
                core_label_y = nicp.core.label['y']
                if nicp.core.is_cyclic_volta:
                    if core_label_x not in dic_x_label:
                        xlabel_set.append(core_label_x)
                        dic_x_label[core_label_x] = 1
                    if core_label_y not in dic_y_label:
                        ylabel_set.append(core_label_y)
                        dic_y_label[core_label_y] = 1
                    if (idx == len(self.multiple_files) - 1):
                        xlabel = ', '.join(xlabel_set)
                        ylabel = ', '.join(ylabel_set)
                else:
                    xlabel = "X ({})".format(core_label_x)
                    ylabel = "Y ({})".format(core_label_y)

            tf.close()
        
        try:
            if global_x_min is not None and global_x_max is not None:
                if any_forward_orientation:
                    plt.xlim(global_x_min, global_x_max)
                else:
                    plt.xlim(global_x_max, global_x_min)
        except Exception:
            pass

        if active_nicp is not None:
            try:
                y_boundary_min, y_boundary_max = active_nicp.plot_overlays(
                    plt, active_y_values, adjust_xlim=False,
                )
                ymin, ymax = plt.gca().get_ylim()
                plt.ylim(min(ymin, y_boundary_min), max(ymax, y_boundary_max))
            except Exception:
                pass

        plt.xlabel(xlabel, fontsize=18)
        plt.ylabel(ylabel, fontsize=18)
        ax = plt.gca()
        if cv_mode:
            ymin, ymax = ax.get_ylim()
            cv_abs_max = max(abs(ymin), abs(ymax))

        if cv_mode and cv_abs_max > 0:
            exp = int(np.floor(np.log10(cv_abs_max))) if cv_abs_max > 0 else 0
            base = (10.0 ** exp) if exp != 0 else 1.0
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _:
                f"{(y / base):.3g}"
            ))
            ax.yaxis.get_offset_text().set_visible(False)
            if exp != 0:
                ax.text(
                    0.0, 1,
                    r"$\times 10^{%d}$" % exp,
                    transform=ax.transAxes,
                    ha='left', va='bottom',
                    fontsize=14,
                    clip_on=False
                )
        plt.grid(False)
        plt.legend()
        tf_img = tempfile.NamedTemporaryFile(suffix='.png')
        plt.savefig(tf_img, format='png')
        tf_img.seek(0)
        plt.clf()
        plt.cla()
        return tf_img
