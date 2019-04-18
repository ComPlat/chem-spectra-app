import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from chem_spectra.model.helper.share import store_str_in_tmp

from chem_spectra.model.converter.jcamp.base import JcampBaseConverter
from chem_spectra.model.converter.jcamp.ni import JcampNIConverter
from chem_spectra.model.converter.jcamp.ms import JcampMSConverter
from chem_spectra.model.converter.ms import MSConverter

from chem_spectra.model.composer.ni import NIComposer
from chem_spectra.model.composer.ms import MSComposer


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
        if not_jcamp:
            return self.ms2composer()

        cv, cp = self.jcamp2cvp()
        return cp


    def ms2composer(self):
        mscv = MSConverter(self.file.core, self.params)
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
