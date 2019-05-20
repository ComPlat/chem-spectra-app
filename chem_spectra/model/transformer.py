import json

from chem_spectra.lib.shared.buffer import store_str_in_tmp
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter
from chem_spectra.lib.converter.jcamp.ms import JcampMSConverter
from chem_spectra.lib.converter.ms import MSConverter
from chem_spectra.lib.composer.ni import NIComposer
from chem_spectra.lib.composer.ms import MSComposer


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
        return cmpsr.tf_jcamp(), cmpsr.tf_img()

    def to_composer(self):
        not_jcamp = self.file.name.split('.')[-1].lower() in ['raw', 'mzml']
        not_jcamp_by_params = self.params['ext'] in ['raw', 'mzml']
        if not_jcamp or not_jcamp_by_params:
            return self.ms2composer()

        cv, cp = self.jcamp2cvp()
        return cp

    def ms2composer(self):
        mscv = MSConverter(self.file, self.params)
        mscp = MSComposer(mscv)
        return mscp

    def jcamp2cvp(self):
        tf = store_str_in_tmp(self.file.core)
        jbcv = JcampBaseConverter(tf.name, self.params)
        tf.close()

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
        r = target.get('output').get('result')
        if not r or not r[0]:
            return False

        tf = store_str_in_tmp(
            json.dumps(target),
            suffix='.json',
        )
        return tf
