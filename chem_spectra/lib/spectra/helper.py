import tempfile
import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from .peaker import SpectraPeaker
from .writer import gen_origin, gen_jcamp_temp


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


def convert_to_img(spPeaker):
    plt.rcParams['figure.figsize'] = [16, 9]
    plt.rcParams['font.size'] = 14
    # PLOT data
    plt.plot(spPeaker.x, spPeaker.y)
    plt.xlim(spPeaker.boundary['x']['max'], spPeaker.boundary['x']['min'])
    # PLOT peaks
    if spPeaker.edit_peaks:
        plt.plot(
            spPeaker.edit_peaks['x'],
            spPeaker.edit_peaks['y'],
            'rd'
        )
    else:
        plt.plot(
            spPeaker.auto_peaks['x'],
            spPeaker.auto_peaks['y'],
            'rd'
        )

    # PLOT label
    plt.xlabel("X ({})".format(spPeaker.label['x']), fontsize=18)
    plt.ylabel("Y ({})".format(spPeaker.label['y']), fontsize=18)
    plt.grid(True)
    # Save
    tf_img = tempfile.NamedTemporaryFile(suffix='.png')
    plt.savefig(tf_img, format='png')
    tf_img.seek(0)
    plt.clf()
    plt.cla()
    return tf_img


def get_jcamp_with_peaks(file, peaks_str=False):
    tf = store_in_tmp(file)
    # identify peaks
    spPeaker = SpectraPeaker(tf.name)
    err = False
    if peaks_str:
        err = spPeaker.write_edit_peak(peaks_str)
    if err or not peaks_str:
        spPeaker.pick_peak()
    # write jcamp
    origin = gen_origin(tf.name)
    tf.close()
    return spPeaker, origin


def convert2jcamp_img(file, peaks_str=False):
    if peaks_str:
        spPeaker, origin = get_jcamp_with_peaks(file, peaks_str)
    else:
        spPeaker, origin = get_jcamp_with_peaks(file)
    tf_jcamp = gen_jcamp_temp(spPeaker, origin)
    tf_img = convert_to_img(spPeaker)
    return tf_jcamp, tf_img


def convert2jcamp(file):
    spPeaker, origin = get_jcamp_with_peaks(file)
    tf_jcamp = gen_jcamp_temp(spPeaker, origin)
    return tf_jcamp


def convert2img(file):
    spPeaker, origin = get_jcamp_with_peaks(file)
    tf_img = convert_to_img(spPeaker)
    return tf_img
