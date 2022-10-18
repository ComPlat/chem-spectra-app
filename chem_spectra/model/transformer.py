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
    def __init__(self, file, molfile=None, params=False):
        self.file = file
        self.molfile = molfile
        self.params = params

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
                    isSimulateNRM = False
                    if self.params and 'simulatenrm' in self.params:
                        isSimulateNRM = self.params['simulatenrm']
                    decorated_jbcv = decorate_sim_property(fbcv, self.molfile, isSimulateNRM)   # noqa: E501
                    # if ((type(decorated_jbcv) is dict) and "invalid_molfile" in decorated_jbcv):
                    #     # return if molfile is invalid
                    #     return None, decorated_jbcv
                    invalid_molfile = False
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
        isSimulateNRM = False
        if params and 'simulatenrm' in params:
            isSimulateNRM = params['simulatenrm']
            
        list_decorated_converters = []
        list_decorated_composers = []
        invalid_molfile = False
        for conv in fid_brucker.data:
            # d_jbcv = decorate_sim_property(conv, self.molfile, isSimulateNRM)   # noqa: E501
            # if ((type(d_jbcv) is dict) and "invalid_molfile" in d_jbcv):
            #     # return if molfile is invalid
            #     return None, d_jbcv

            decorated_jbcv = decorate_sim_property(conv, self.molfile, isSimulateNRM)   # noqa: E501
            
            if ((type(decorated_jbcv) is dict) and "invalid_molfile" in decorated_jbcv):
                invalid_molfile = True
                final_decorated_jbcv = decorated_jbcv['origin_jbcv']
            else:
                final_decorated_jbcv = decorated_jbcv
            
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
            # isSimulateNRM = self.params['simulatenrm']
            isSimulateNRM = False
            if self.params and 'simulatenrm' in self.params:
                isSimulateNRM = self.params['simulatenrm']
            # d_jbcv = decorate_sim_property(jbcv, self.molfile, isSimulateNRM)
            # if ((type(d_jbcv) is dict) and "invalid_molfile" in d_jbcv):
            #     # return if molfile is invalid
            #     return None, d_jbcv
            decorated_jbcv = decorate_sim_property(jbcv, self.molfile, isSimulateNRM)   # noqa: E501
            
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
        nicp = NIComposer(converter)
        tf_jcamp = nicp.tf_jcamp()
        return tf_jcamp
