import pytest
import io

from werkzeug.datastructures import FileStorage
from chem_spectra.lib.spectra.helper import (
    convert_jcamp_temp, create_sp_carrier
)
from chem_spectra.lib.spectra.writer import (
    build_block_meta
)

from chem_spectra.lib.spectra.carrier import SpectraCarrier


target_dir = './tests/fixtures/'
source_dir = 'source/'
result_dir = 'result/'

IR_dx = 'IR.dx'
H1_dx = '1H.dx'
C13_CPD_dx = '13C-CPD.dx'
C13_DEPT135_dx = '13C-DEPT135.dx'
SVS_790A_13C_jdx = 'SVS-790A_13C.jdx'


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# generate datatable
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def __fixture_path(orig_filename):
    return target_dir + source_dir + orig_filename


def __generated_jcamp_temp(path, params=False):
    with open(path, 'rb') as f:
        file = FileStorage(f)
        sp_carrier = create_sp_carrier(file, params)
        tmp = convert_jcamp_temp(sp_carrier)
    return sp_carrier, tmp


def __target_peaks_meta(filename):
    meta_target = open(target_dir + result_dir + filename, 'rb')
    return meta_target.read().decode('utf-8')


def __tolerable(one, two, ref, defined_tolerance=None):
    tolerance = defined_tolerance or 1 / 100000
    err = abs((one - two) / ref)
    return err < tolerance


def __is_match(ori, nxt, ref, count_target, tolerance=None):
    count_actual = 0
    for i in range(count_target):
        if ori[i] == nxt[i]:
            count_actual += 1
        elif __tolerable(ori[i], nxt[i], ref, tolerance):
            count_actual += 1
    return count_target == count_actual


def test_datatable_IR():
    sp_car_ori, ori_tmp = __generated_jcamp_temp(
        __fixture_path(IR_dx)
    )
    sp_car_nxt, nxt_tmp = __generated_jcamp_temp(ori_tmp.name)
    sp_car_ano, ano_tmp = __generated_jcamp_temp(nxt_tmp.name)
    total_count = sp_car_ori.y.shape[0]
    ori_bd = sp_car_ori.boundary
    ref = abs(ori_bd['x']['max'] - ori_bd['x']['min'])

    assert __is_match(sp_car_ori.x, sp_car_nxt.x, ref, total_count)
    assert __is_match(sp_car_ori.y, sp_car_nxt.y, ref, total_count)
    assert __is_match(sp_car_ori.x, sp_car_ano.x, ref, total_count)
    assert __is_match(sp_car_ori.y, sp_car_ano.y, ref, total_count)


def test_datatable_1H():
    sp_car_ori, ori_tmp = __generated_jcamp_temp(
        __fixture_path(H1_dx)
    )
    sp_car_nxt, nxt_tmp = __generated_jcamp_temp(ori_tmp.name)
    sp_car_ano, ano_tmp = __generated_jcamp_temp(nxt_tmp.name)
    total_count = sp_car_ori.y.shape[0]
    ori_bd = sp_car_ori.boundary
    ref = abs(ori_bd['x']['max'] - ori_bd['x']['min'])

    assert __is_match(sp_car_ori.x, sp_car_nxt.x, ref, total_count)
    assert __is_match(sp_car_ori.y, sp_car_nxt.y, ref, total_count)
    assert __is_match(sp_car_ori.x, sp_car_ano.x, ref, total_count)
    assert __is_match(sp_car_ori.y, sp_car_ano.y, ref, total_count)


def test_datatable_13C_CPD_dx():
    sp_car_ori, ori_tmp = __generated_jcamp_temp(
        __fixture_path(C13_CPD_dx)
    )
    sp_car_nxt, nxt_tmp = __generated_jcamp_temp(ori_tmp.name)
    sp_car_ano, ano_tmp = __generated_jcamp_temp(nxt_tmp.name)
    total_count = sp_car_ori.y.shape[0]
    ori_bd = sp_car_ori.boundary
    ref = abs(ori_bd['x']['max'] - ori_bd['x']['min'])

    assert __is_match(sp_car_ori.x, sp_car_nxt.x, ref, total_count)
    assert __is_match(sp_car_ori.y, sp_car_nxt.y, ref, total_count)
    assert __is_match(sp_car_ori.x, sp_car_ano.x, ref, total_count)
    assert __is_match(sp_car_ori.y, sp_car_ano.y, ref, total_count)


def test_datatable_13C_DEPT135():
    sp_car_ori, ori_tmp = __generated_jcamp_temp(
        __fixture_path(C13_DEPT135_dx)
    )
    sp_car_nxt, nxt_tmp = __generated_jcamp_temp(ori_tmp.name)
    sp_car_ano, ano_tmp = __generated_jcamp_temp(nxt_tmp.name)
    total_count = sp_car_ori.y.shape[0]
    ori_bd = sp_car_ori.boundary
    ref = abs(ori_bd['x']['max'] - ori_bd['x']['min'])

    assert __is_match(sp_car_ori.x, sp_car_nxt.x, ref, total_count)
    assert __is_match(sp_car_ori.y, sp_car_nxt.y, ref, total_count)
    assert __is_match(sp_car_ori.x, sp_car_ano.x, ref, total_count)
    assert __is_match(sp_car_ori.y, sp_car_ano.y, ref, total_count)


def test_datatable_SVS_790A_13C_jdx():
    sp_car_ori, ori_tmp = __generated_jcamp_temp(
        __fixture_path(SVS_790A_13C_jdx)
    )
    sp_car_nxt, nxt_tmp = __generated_jcamp_temp(ori_tmp.name)
    sp_car_ano, ano_tmp = __generated_jcamp_temp(nxt_tmp.name)
    total_count = sp_car_ori.y.shape[0]
    ori_bd = sp_car_ori.boundary
    ref = abs(ori_bd['x']['max'] - ori_bd['x']['min'])

    assert __is_match(sp_car_ori.x, sp_car_nxt.x, ref, total_count)
    assert __is_match(sp_car_ori.y, sp_car_nxt.y, ref, total_count, 1/10000)
    assert __is_match(sp_car_ori.x, sp_car_ano.x, ref, total_count)
    assert __is_match(sp_car_ori.y, sp_car_ano.y, ref, total_count, 1/10000)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# compare full file
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def __is_same(ones, twos):
    inconsistent = 0
    for i in range(len(ones)):
        if ones[i] != twos[i] and ones[i][:7] != '##TITLE':
            inconsistent += 1
    return inconsistent == 0


def test_compare_datatable_13C_DEPT135():
    sp_car_ori, ori_tmp = __generated_jcamp_temp(
        __fixture_path(C13_DEPT135_dx)
    )
    sp_car_nxt, nxt_tmp = __generated_jcamp_temp(ori_tmp.name)
    sp_car_ano, ano_tmp = __generated_jcamp_temp(nxt_tmp.name)

    meta_ori = build_block_meta(sp_car_ori)
    meta_nxt = build_block_meta(sp_car_nxt)
    meta_ano = build_block_meta(sp_car_ano)

    assert(__is_same(meta_ori, meta_nxt))
    assert(__is_same(meta_ori, meta_ano))
