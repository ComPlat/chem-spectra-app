import pytest
import io

from werkzeug.datastructures import FileStorage
from chem_spectra.lib.spectra.helper import (
    convert2jcamp_img, create_sp_carrier
)
from chem_spectra.lib.spectra.writer import (
    build_block_meta
)


target_dir = './tests/fixtures/'
source_dir = 'source/'
result_dir = 'result/'

IR_dx = 'IR.dx'
H1_dx = '1H.dx'
C13_CPD_dx = '13C-CPD.dx'
C13_DEPT135_dx = '13C-DEPT135.dx'
SVS_790A_13C_jdx = 'SVS-790A_13C.jdx'

meta_IR_dx = 'meta_IR'
meta_H1_dx = 'meta_1H'
meta_C13_CPD_dx = 'meta_13C-CPD'
meta_C13_DEPT135_dx = 'meta_13C-DEPT135'
meta_SVS_790A_13C_jdx = 'meta_SVS-790A_13C'

separator = '$$ === CHEMSPECTRA PEAK ASSIGNMENTS AUTO ==='


def __generated_peaks_meta(orig_filename, params=False):
    with open(target_dir + source_dir + orig_filename, 'rb') as f:
        file = FileStorage(f)
        sp_carrier = create_sp_carrier(file, params)
        meta_jcamp = build_block_meta(sp_carrier)

    parts = ''.join(meta_jcamp).split(separator)
    assert len(parts) == 2

    meta_content = parts[1]
    return meta_content


def __target_peaks_meta(filename):
    meta_target = open(target_dir + result_dir + filename, 'rb')
    return meta_target.read().decode('utf-8')


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# generate peaks + origin files
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def test_meta_IR():
    meta_content = __generated_peaks_meta(IR_dx)
    meta_target = __target_peaks_meta(meta_IR_dx)
    assert meta_content == meta_target


def test_meta_1H():
    meta_content = __generated_peaks_meta(H1_dx)
    meta_target = __target_peaks_meta(meta_H1_dx)
    assert meta_content == meta_target


def test_meta_13C_CPD_dx():
    meta_content = __generated_peaks_meta(C13_CPD_dx)
    meta_target = __target_peaks_meta(meta_C13_CPD_dx)
    assert meta_content == meta_target


def test_meta_13C_DEPT135():
    meta_content = __generated_peaks_meta(C13_DEPT135_dx)
    meta_target = __target_peaks_meta(meta_C13_DEPT135_dx)
    assert meta_content == meta_target


def test_meta_SVS_790A_13C_jdx():
    meta_content = __generated_peaks_meta(SVS_790A_13C_jdx)
    meta_target = __target_peaks_meta(meta_SVS_790A_13C_jdx)
    assert meta_content == meta_target


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# generate peaks + edit files
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def test_edit_meta_IR():
    meta_content = __generated_peaks_meta('edit/edit_' + IR_dx)
    meta_target = __target_peaks_meta('edit/edit_' + meta_IR_dx)
    assert meta_content == meta_target


def test_edit_meta_1H():
    meta_content = __generated_peaks_meta('edit/edit_' + H1_dx)
    meta_target = __target_peaks_meta('edit/edit_' + meta_H1_dx)
    assert meta_content == meta_target


def test_edit_meta_13C_DEPT135():
    meta_content = __generated_peaks_meta('edit/edit_' + C13_DEPT135_dx)
    meta_target = __target_peaks_meta('edit/edit_' + meta_C13_DEPT135_dx)
    assert meta_content == meta_target


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# edit peaks + origin files
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
params_1 = {
    'peaks_str': '1.1,1.112#2.2,2.224#3.3,3.336',
}

def test_params_meta_IR_dx():
    ps_meta_content = __generated_peaks_meta(IR_dx, params_1)
    meta_target = __target_peaks_meta('ps/ps_' + meta_IR_dx)
    assert ps_meta_content == meta_target


def test_params_meta_1H():
    ps_meta_content = __generated_peaks_meta(H1_dx, params_1)
    meta_target = __target_peaks_meta('ps/ps_' + meta_H1_dx)
    assert ps_meta_content == meta_target


def test_params_meta_13C_DEPT135():
    ps_meta_content = __generated_peaks_meta(C13_DEPT135_dx, params_1)
    meta_target = __target_peaks_meta('ps/ps_' + meta_C13_DEPT135_dx)
    assert ps_meta_content == meta_target


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# edit peaks + auto files
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def test_auto_params_meta_IR_dx():
    auto_meta_content = __generated_peaks_meta('auto/auto_' + IR_dx, params_1)
    meta_target = __target_peaks_meta('auto/auto_ps_' + meta_IR_dx)
    assert auto_meta_content == meta_target


def test_auto_params_meta_1H():
    auto_meta_content = __generated_peaks_meta('auto/auto_' + H1_dx, params_1)
    meta_target = __target_peaks_meta('auto/auto_ps_' + meta_H1_dx)
    assert auto_meta_content == meta_target


def test_auto_params_auto_meta_13C_DEPT135():
    auto_meta_content = __generated_peaks_meta('auto/auto_' + C13_DEPT135_dx, params_1)
    meta_target = __target_peaks_meta('auto/auto_ps_' + meta_C13_DEPT135_dx)
    assert auto_meta_content == meta_target


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# edit peaks + edit files
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def test_edit_params_meta_IR_dx():
    edit_meta_content = __generated_peaks_meta('edit/edit_' + IR_dx, params_1)
    meta_target = __target_peaks_meta('auto/auto_ps_' + meta_IR_dx)
    assert edit_meta_content == meta_target


def test_auto_params_meta_1H():
    auto_meta_content = __generated_peaks_meta('edit/edit_' + H1_dx, params_1)
    meta_target = __target_peaks_meta('auto/auto_ps_' + meta_H1_dx)
    assert auto_meta_content == meta_target


def test_auto_params_auto_meta_13C_DEPT135():
    auto_meta_content = __generated_peaks_meta('edit/edit_' + C13_DEPT135_dx, params_1)
    meta_target = __target_peaks_meta('auto/auto_ps_' + meta_C13_DEPT135_dx)
    assert auto_meta_content == meta_target


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# given EMPTY peaks string
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
empty_params = {
    'peaks_str': '',
}


def test_empty_params_meta_IR_dx():
    empty_meta_content = __generated_peaks_meta(IR_dx, empty_params)
    meta_target = __target_peaks_meta(meta_IR_dx)
    assert empty_meta_content == meta_target


def test_empty_params_meta_1H():
    empty_meta_content = __generated_peaks_meta(H1_dx, empty_params)
    meta_target = __target_peaks_meta(meta_H1_dx)
    assert empty_meta_content == meta_target


def test_empty_params_meta_13C_DEPT135():
    empty_meta_content = __generated_peaks_meta(C13_DEPT135_dx, empty_params)
    meta_target = __target_peaks_meta(meta_C13_DEPT135_dx)
    assert empty_meta_content == meta_target


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# edit peaks / shift + edit files
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
params_2 = {
    'peaks_str': '1.1,1.112#2.2,2.224#3.3,3.336',
    'select_x': '40.0',
    'ref_name': 'Actic acid-d4 (s)',
    'ref_value': '160.0',
}


def test_shift_edit_params_meta_1H():
    auto_meta_content = __generated_peaks_meta('edit/edit_' + H1_dx, params_2)
    meta_target = __target_peaks_meta('shift/shift_edit_' + meta_H1_dx)
    assert auto_meta_content == meta_target


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# inherit sample description
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
params_3 = {
}


def test_inherit_wihtout_params_1H():
    inherit_meta_content = __generated_peaks_meta('edit/inherit_' + H1_dx, params_3)
    meta_target = __target_peaks_meta('edit/inherit_' + meta_H1_dx)
    assert inherit_meta_content == meta_target
