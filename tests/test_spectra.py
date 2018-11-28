import pytest
import io

from werkzeug.datastructures import FileStorage
from chem_spectra.lib.spectra.helper import (
    convert2jcamp_img, convert2jcamp, convert2img,
    get_jcamp_with_peaks
)
from chem_spectra.lib.spectra.writer import (
    gen_jcamp_meta
)


target_dir = './tests/fixtures/'
file_jdx = '13C-DEPT135.dx'
another_jdx = 'IR.dx'
result_jdx = 'result.jdx'
result_png = 'result.png'
result_meta = 'result_meta.jdx'
separator = '$$ === CHEMSPECTRA PEAK ASSIGNMENTS ALL ==='


def test_meta_jcamp():
    with open(target_dir + file_jdx, 'rb') as f:
        file = FileStorage(f)
        spPeaker, origin = get_jcamp_with_peaks(file)
        meta_jcamp = gen_jcamp_meta(spPeaker, origin)

    meta_content = ''.join(meta_jcamp).split(separator)[1]
    meta_target = open(target_dir + result_meta, 'rb')
    assert meta_content == meta_target.read().decode('utf-8')


# def test_convert2jcamp():
#     with open(target_dir + file_jdx, 'rb') as f:
#         file = FileStorage(f)
#         tf_jcamp = convert2jcamp(file)

#     f_jdx = open(target_dir + result_jdx, 'rb')

#     assert tf_jcamp.read() == f_jdx.read()


# def test_convert2img():
#     with open(target_dir + file_jdx, 'rb') as f:
#         file = FileStorage(f)
#         tf_img = convert2img(file)

#     f_png = open(target_dir + result_png, 'rb')

#     assert tf_img.read() == f_png.read()


# def test_convert2jcamp_img():
#     with open(target_dir + file_jdx, 'rb') as f:
#         file = FileStorage(f)
#         tf_jcamp, tf_img = convert2jcamp_img(file)

#     f_jdx = open(target_dir + result_jdx, 'rb')
#     f_png = open(target_dir + result_png, 'rb')

#     assert tf_jcamp.read() == f_jdx.read()
#     assert tf_img.read() == f_png.read()

#     # - - - - -
#     with open(target_dir + another_jdx, 'rb') as f:
#         file = FileStorage(f)
#         tf_jcamp, tf_img = convert2jcamp_img(file)

#     f_jdx = open(target_dir + result_jdx, 'rb')
#     f_png = open(target_dir + result_png, 'rb')

#     assert tf_jcamp.read() != f_jdx.read()
#     assert tf_img.read() != f_png.read()
