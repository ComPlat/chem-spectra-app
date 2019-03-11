import tempfile
import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from chem_spectra.lib.spectra.carrier import SpectraCarrier
from chem_spectra.lib.spectra.writer import convert_jcamp_temp
from chem_spectra.lib.ms.raw_converter import RawConverter
from chem_spectra.lib.ms.ms_composer import MsComposer

ALLOWED_EXTENSIONS = set(['dx', 'jdx'])


def allowed_file(file):
    ext = file.filename.split('.')[-1].lower()
    return ext in ALLOWED_EXTENSIONS


def store_in_tmp(file):
    byteContent = file.stream.read()
    tf = tempfile.NamedTemporaryFile()
    with open(tf.name, 'w') as f:
        tf.write(byteContent)
    return tf


def convert_to_img(sp_carrier):
    plt.rcParams['figure.figsize'] = [16, 9]
    plt.rcParams['font.size'] = 14
    # PLOT data
    plt.plot(sp_carrier.x, sp_carrier.y)
    plt.xlim(sp_carrier.boundary['x']['max'], sp_carrier.boundary['x']['min'])
    # PLOT peaks
    if sp_carrier.edit_peaks:
        plt.plot(
            sp_carrier.edit_peaks['x'],
            sp_carrier.edit_peaks['y'],
            'rd'
        )
    elif sp_carrier.auto_peaks:
        plt.plot(
            sp_carrier.auto_peaks['x'],
            sp_carrier.auto_peaks['y'],
            'rd'
        )

    # PLOT label
    plt.xlabel("X ({})".format(sp_carrier.label['x']), fontsize=18)
    plt.ylabel("Y ({})".format(sp_carrier.label['y']), fontsize=18)
    plt.grid(False)
    # Save
    tf_img = tempfile.NamedTemporaryFile(suffix='.png')
    plt.savefig(tf_img, format='png')
    tf_img.seek(0)
    plt.clf()
    plt.cla()
    return tf_img


def create_sp_carrier(file, params=False):
    tf = store_in_tmp(file)
    # identify peaks, etc
    sp_carrier = SpectraCarrier(tf.name, params)
    tf.close()
    return sp_carrier


def convert2jcamp_img(file, params=False):
    sp_carrier = create_sp_carrier(file, params)
    tf_jcamp = convert_jcamp_temp(sp_carrier)
    tf_img = convert_to_img(sp_carrier)
    return tf_jcamp, tf_img


def convert2jcamp(file, params=False):
    sp_carrier = create_sp_carrier(file, params)
    tf_jcamp = convert_jcamp_temp(sp_carrier)
    return tf_jcamp


def convert2img(file, params=False):
    sp_carrier = create_sp_carrier(file, params)
    tf_img = convert_to_img(sp_carrier)
    return tf_img


def convertRaw2jcamp_img(file, exact_mz=0):
    rc = RawConverter(file, exact_mz)
    mc = MsComposer(rc)
    return mc.tf_jcamp(), mc.tf_img()
