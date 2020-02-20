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


def search_fid(td):
    try:
        for one in os.listdir(td):
            targets = glob.glob('{}/{}/fid'.format(td, one))
            if len(targets) > 0:
                return '{}/{}'.format(td, one)
            if os.path.isdir('{}/{}'.format(td, one)):
                for two in os.listdir('{}/{}'.format(td, one)):
                    targets = glob.glob('{}/{}/{}/fid'.format(td, one, two))[0]
                    if len(targets) > 0:
                        return '{}/{}/{}'.format(td, one, two)
        return False
    except:
        return False


class TransformerModel:
    def __init__(self, file, params=False):
        self.file = file
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
        is_raw_mzml = self.file.name.split('.')[-1].lower() in ['raw', 'mzml']
        is_cdf = self.file.name.split('.')[-1].lower() in ['cdf']
        is_zip = self.file.name.split('.')[-1].lower() in ['zip']
        is_raw_mzml_by_params = self.params['ext'] in ['raw', 'mzml']
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
            target_dir = search_fid(td)
            if not target_dir:
                return False, False
            fbcv = FidBaseConverter(target_dir, self.params, self.file.name)
            if not fbcv:
                return False, False
        # assume NMR only
        nicv = JcampNIConverter(fbcv)
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
            nicv = JcampNIConverter(jbcv)
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
