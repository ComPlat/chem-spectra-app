import pytest
import io

from werkzeug.datastructures import FileStorage
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.model.transformer import TransformerModel as TraModel

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
        file = FileContainer(FileStorage(f))
        nicv, nicp = TraModel(file, params).jcamp2cvp()
        jcamp = nicp.tf_jcamp()
    return nicv, nicp, jcamp


def __target_peaks_meta(filename):
    meta_target = open(target_dir + result_dir + filename, 'rb')
    return meta_target.read().decode('utf-8', errors='ignore')


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
    nicv_ori, nicp_ori, jcamp_ori = __generated_jcamp_temp(
        __fixture_path(IR_dx)
    )
    nicv_nxt, nicp_nxt, jcamp_nxt = __generated_jcamp_temp(jcamp_ori.name)
    nicv_ano, nicp_ano, jcamp_ano = __generated_jcamp_temp(jcamp_nxt.name)
    total_count = nicv_ori.ys.shape[0]
    ori_bd = nicv_ori.boundary
    ref = abs(ori_bd['x']['max'] - ori_bd['x']['min'])

    assert __is_match(nicv_ori.xs, nicv_nxt.xs, ref, total_count)
    assert __is_match(nicv_ori.ys, nicv_nxt.ys, ref, total_count)
    assert __is_match(nicv_ori.xs, nicv_nxt.xs, ref, total_count)
    assert __is_match(nicv_ori.ys, nicv_nxt.ys, ref, total_count)


def test_datatable_1H():
    nicv_ori, nicp_ori, jcamp_ori = __generated_jcamp_temp(
        __fixture_path(H1_dx)
    )
    nicv_nxt, nicp_nxt, jcamp_nxt = __generated_jcamp_temp(jcamp_ori.name)
    nicv_ano, nicp_ano, jcamp_ano = __generated_jcamp_temp(jcamp_nxt.name)
    total_count = nicv_ori.ys.shape[0]
    ori_bd = nicv_ori.boundary
    ref = abs(ori_bd['x']['max'] - ori_bd['x']['min'])

    assert __is_match(nicv_ori.xs, nicv_nxt.xs, ref, total_count)
    assert __is_match(nicv_ori.ys, nicv_nxt.ys, ref, total_count)
    assert __is_match(nicv_ori.xs, nicv_ano.xs, ref, total_count)
    assert __is_match(nicv_ori.ys, nicv_ano.ys, ref, total_count)


def test_datatable_13C_CPD_dx():
    nicv_ori, nicp_ori, jcamp_ori = __generated_jcamp_temp(
        __fixture_path(C13_CPD_dx)
    )
    nicv_nxt, nicp_nxt, jcamp_nxt = __generated_jcamp_temp(jcamp_ori.name)
    nicv_ano, nicp_ano, jcamp_ano = __generated_jcamp_temp(jcamp_nxt.name)
    total_count = nicv_ori.ys.shape[0]
    ori_bd = nicv_ori.boundary
    ref = abs(ori_bd['x']['max'] - ori_bd['x']['min'])

    assert __is_match(nicv_ori.xs, nicv_nxt.xs, ref, total_count)
    assert __is_match(nicv_ori.ys, nicv_nxt.ys, ref, total_count)
    assert __is_match(nicv_ori.xs, nicv_ano.xs, ref, total_count)
    assert __is_match(nicv_ori.ys, nicv_ano.ys, ref, total_count)


def test_datatable_13C_DEPT135():
    nicv_ori, nicp_ori, jcamp_ori = __generated_jcamp_temp(
        __fixture_path(C13_DEPT135_dx)
    )
    nicv_nxt, nicp_nxt, jcamp_nxt = __generated_jcamp_temp(jcamp_ori.name)
    nicv_ano, nicp_ano, jcamp_ano = __generated_jcamp_temp(jcamp_nxt.name)
    total_count = nicv_ori.ys.shape[0]
    ori_bd = nicv_ori.boundary
    ref = abs(ori_bd['x']['max'] - ori_bd['x']['min'])

    assert __is_match(nicv_ori.xs, nicv_nxt.xs, ref, total_count)
    assert __is_match(nicv_ori.ys, nicv_nxt.ys, ref, total_count)
    assert __is_match(nicv_ori.xs, nicv_ano.xs, ref, total_count)
    assert __is_match(nicv_ori.ys, nicv_ano.ys, ref, total_count)


def test_datatable_SVS_790A_13C_jdx():
    nicv_ori, nicp_ori, jcamp_ori = __generated_jcamp_temp(
        __fixture_path(SVS_790A_13C_jdx)
    )
    nicv_nxt, nicp_nxt, jcamp_nxt = __generated_jcamp_temp(jcamp_ori.name)
    nicv_ano, nicp_ano, jcamp_ano = __generated_jcamp_temp(jcamp_nxt.name)
    total_count = nicv_ori.ys.shape[0]
    ori_bd = nicv_ori.boundary
    ref = abs(ori_bd['x']['max'] - ori_bd['x']['min'])

    assert __is_match(nicv_ori.xs, nicv_nxt.xs, ref, total_count)
    assert __is_match(nicv_ori.ys, nicv_nxt.ys, ref, total_count, 1/10000)
    assert __is_match(nicv_ori.xs, nicv_ano.xs, ref, total_count)
    assert __is_match(nicv_ori.ys, nicv_ano.ys, ref, total_count, 1/10000)


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
    nicv_ori, nicp_ori, jcamp_ori = __generated_jcamp_temp(
        __fixture_path(C13_DEPT135_dx)
    )
    nicv_nxt, nicp_nxt, jcamp_nxt = __generated_jcamp_temp(jcamp_ori.name)
    nicv_ano, nicp_ano, jcamp_ano = __generated_jcamp_temp(jcamp_nxt.name)

    assert(__is_same(nicp_nxt.meta, nicp_ano.meta))
