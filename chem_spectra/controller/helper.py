import tempfile
import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from chem_spectra.model.converter.jcamp.base import JcampBaseConverter
from chem_spectra.model.converter.jcamp.ni import JcampNIConverter
from chem_spectra.model.converter.jcamp.ms import JcampMSConverter
from chem_spectra.model.converter.ms import MSConverter

from chem_spectra.model.composer.ni import NIComposer
from chem_spectra.model.composer.ms import MSComposer


ALLOWED_EXTENSIONS = set(['dx', 'jdx', 'raw', 'mzml'])


def allowed_file(file):
    ext = file.filename.split('.')[-1].lower()
    return ext in ALLOWED_EXTENSIONS


def store_in_tmp(file):
    byteContent = file.stream.read()
    tf = tempfile.NamedTemporaryFile()
    with open(tf.name, 'w') as f:
        tf.write(byteContent)
    return tf


def convert2jcamp(file, params=False):
    cmpsr = to_composer(file, params)
    return cmpsr.tf_jcamp()


def convert2img(file, params=False):
    cmpsr = to_composer(file, params)
    return cmpsr.tf_img()


def convert2jcamp_img(file, params=False):
    cmpsr = to_composer(file, params)
    return cmpsr.tf_jcamp(), cmpsr.tf_img()


def to_composer(file, params):
    not_jcamp = file.filename.split('.')[-1].lower() in ['raw', 'mzml']
    if not_jcamp:
        return ms2composer(file, params)

    cv, cp = jcamp2cvp(file, params)
    return cp


def ms2composer(file, params=False):
    mscv = MSConverter(file, params)
    mscp = MSComposer(mscv)
    return mscp


def jcamp2cvp(file, params=False):
    tf = store_in_tmp(file)
    jbcv = JcampBaseConverter(tf.name, params)
    tf.close()

    if jbcv.typ == 'MS':
        mscv = JcampMSConverter(jbcv)
        mscp = MSComposer(mscv)
        return mscv, mscp
    else:
        nicv = JcampNIConverter(jbcv)
        nicp = NIComposer(nicv)
        return nicv, nicp
