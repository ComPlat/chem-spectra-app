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
        if not target_dir:
            target_dir = find_dir(td, 'ser')
        return target_dir, has_processed_files
    except:     # noqa: E722
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
        if converter.is_2d == True:
            return None
        nicp = NIComposer(converter)
        tf_jcamp = nicp.tf_jcamp()
        return tf_jcamp
      
    def tf_combine(self, list_file_names=None):
        if not self.multiple_files:
            return False
          
        plt.rcParams['figure.figsize'] = [16, 9]
        plt.rcParams['font.size'] = 14
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
        
        plt.xlabel(xlabel, fontsize=18)
        plt.ylabel(ylabel, fontsize=18)
        plt.legend()
        tf_img = tempfile.NamedTemporaryFile(suffix='.png')
        plt.savefig(tf_img, format='png')
        tf_img.seek(0)
        plt.clf()
        plt.cla()
        return tf_img
