import json
import zipfile
import tempfile
import glob
import os

from chem_spectra.lib.shared.buffer import store_str_in_tmp, store_byte_in_tmp
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter
from chem_spectra.lib.converter.jcamp.ms import JcampMSConverter
from chem_spectra.lib.converter.cdf.base import CdfBaseConverter
from chem_spectra.lib.converter.cdf.ms import CdfMSConverter
from chem_spectra.lib.converter.fid.base import FidBaseConverter
from chem_spectra.lib.converter.ms import MSConverter
from chem_spectra.lib.composer.ni import NIComposer
from chem_spectra.lib.composer.ms import MSComposer

from chem_spectra.model.concern.property import decorate_sim_property


def find_dir(path, name):
    for root, _, files in os.walk(path):
        if name in files:
            return os.path.join(root)
    return False


def search_brucker_binary(td):
    try:
        target_dir = find_dir(td, 'fid')
        if not target_dir:
            target_dir = find_dir(td, 'ser')
        return target_dir
    except:
        return False


class TransformerModel:
    def __init__(self, file, molfile=None, params=False):
        self.file = file
        self.molfile  = molfile
        self.params = params

    def convert2jcamp(self):
        cmpsr = self.to_composer()
        return cmpsr.tf_jcamp()

    def convert2img(self):
        cmpsr = self.to_composer()
        return cmpsr.tf_img()

    def convert2jcamp_img(self):
        cmpsr = self.to_composer()
        if not cmpsr:
            return False, False
        return cmpsr.tf_jcamp(), cmpsr.tf_img()

    def to_composer(self):
        is_raw_mzml = self.file.name.split('.')[-1].lower() in ['raw', 'mzml', 'mzxml']
        is_cdf = self.file.name.split('.')[-1].lower() in ['cdf']
        is_zip = self.file.name.split('.')[-1].lower() in ['zip']
        is_raw_mzml_by_params = self.params['ext'] in ['raw', 'mzml', 'mzxml']
        is_cdf_by_params = self.params['ext'] in ['cdf']
        is_zip_by_params = self.params['ext'] in ['zip']
        if is_raw_mzml or is_raw_mzml_by_params:
            return self.ms2composer()
        if is_cdf or is_cdf_by_params:
            _, cp = self.cdf2cvp()
            return cp
        if is_zip or is_zip_by_params:
            _, cp = self.zip2cvp()
            return cp
        else:
            _, cp = self.jcamp2cvp()
            return cp

    def to_converter(self):
        is_raw_mzml = self.file.name.split('.')[-1].lower() in ['raw', 'mzml', 'mzxml']
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
            target_dir = search_brucker_binary(td)
            if not target_dir:
                return False, False
            fbcv = FidBaseConverter(target_dir, self.params, self.file.name)
            if not fbcv:
                return False, False
        # assume NMR only
        isSimulateNRM = self.params['simulatenrm']
        d_jbcv = decorate_sim_property(fbcv, self.molfile, isSimulateNRM)
        if ((type(d_jbcv) is dict) and "invalid_molfile" in d_jbcv):
            #return if molfile is invalid
            return None, d_jbcv
        nicv = JcampNIConverter(d_jbcv)
        nicp = NIComposer(nicv)
        return nicv, nicp

    def jcamp2cvp(self):
        tf = store_str_in_tmp(self.file.core)
        jbcv = JcampBaseConverter(tf.name, self.params)
        tf.close()
        # conversion
        if jbcv.typ == 'MS':
            mscv = JcampMSConverter(jbcv)
            mscp = MSComposer(mscv)
            return mscv, mscp
        else:
            isSimulateNRM = self.params['simulatenrm']
            d_jbcv = decorate_sim_property(jbcv, self.molfile, isSimulateNRM)
            if ((type(d_jbcv) is dict) and "invalid_molfile" in d_jbcv):
                #return if molfile is invalid
                return None, d_jbcv
            nicv = JcampNIConverter(d_jbcv)
            nicp = NIComposer(nicv)
            return nicv, nicp

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
