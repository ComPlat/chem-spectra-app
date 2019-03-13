import tempfile
import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from chem_spectra.model.converter.nmr_ir import NmrIrConverter
from chem_spectra.model.converter.raw import RawConverter
from chem_spectra.model.composer.nmr_ir import NmrIrComposer
from chem_spectra.model.composer.ms import MsComposer


ALLOWED_EXTENSIONS = set(['dx', 'jdx', 'raw'])


def allowed_file(file):
    ext = file.filename.split('.')[-1].lower()
    return ext in ALLOWED_EXTENSIONS


def store_in_tmp(file):
    byteContent = file.stream.read()
    tf = tempfile.NamedTemporaryFile()
    with open(tf.name, 'w') as f:
        tf.write(byteContent)
    return tf


def create_nicv(file, params=False):
    tf = store_in_tmp(file)
    nicv = NmrIrConverter(tf.name, params)
    tf.close()
    return nicv


def convert2jcamp(file, params=False):
    is_raw = file.filename.split('.')[-1].lower() == 'raw'
    if is_raw:
        return raw2jcamp_img(file, params)[0]

    nicv = create_nicv(file, params)
    nicp = NmrIrComposer(nicv)
    return nicp.tf_jcamp()


def convert2img(file, params=False):
    is_raw = file.filename.split('.')[-1].lower() == 'raw'
    if is_raw:
        return raw2jcamp_img(file, params)[1]

    nicv = create_nicv(file, params)
    nicp = NmrIrComposer(nicv)
    return nicp.tf_img()


def convert2jcamp_img(file, params=False):
    is_raw = file.filename.split('.')[-1].lower() == 'raw'
    if is_raw:
        return raw2jcamp_img(file, params)

    nicv = create_nicv(file, params)
    nicp = NmrIrComposer(nicv)
    return nicp.tf_jcamp(), nicp.tf_img()


def raw2jcamp_img(file, params=False):
    rc = RawConverter(file, params)
    mc = MsComposer(rc)
    return mc.tf_jcamp(), mc.tf_img()
