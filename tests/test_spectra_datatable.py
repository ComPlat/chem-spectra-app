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
JPK_948_jdx = 'JPK-948.jdx'
MS_dx = '/ms/ms_v6.dx'


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
        molfile = FileContainer(FileStorage(None))
        nicv, nicp, _ = TraModel(file, molfile, params).jcamp2cvp()
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
    total_count = nicv_ori.ys.shape[0]
    ori_bd = nicv_ori.boundary
    ref = abs(ori_bd['x']['max'] - ori_bd['x']['min'])

    assert __is_match(nicv_ori.xs, nicv_nxt.xs, ref, total_count)
    assert __is_match(nicv_ori.ys, nicv_nxt.ys, ref, total_count)


def test_datatable_1H():
    nicv_ori, nicp_ori, jcamp_ori = __generated_jcamp_temp(
        __fixture_path(H1_dx)
    )
    nicv_nxt, nicp_nxt, jcamp_nxt = __generated_jcamp_temp(jcamp_ori.name)
    total_count = nicv_ori.ys.shape[0]
    ori_bd = nicv_ori.boundary
    ref = 1

    assert __is_match(nicv_ori.xs, nicv_nxt.xs, ref, total_count)
    assert __is_match(nicv_ori.ys, nicv_nxt.ys, ref, total_count)


def test_datatable_13C_CPD_dx():
    nicv_ori, nicp_ori, jcamp_ori = __generated_jcamp_temp(
        __fixture_path(C13_CPD_dx)
    )
    nicv_nxt, nicp_nxt, jcamp_nxt = __generated_jcamp_temp(jcamp_ori.name)
    total_count = nicv_ori.ys.shape[0]
    ori_bd = nicv_ori.boundary
    ref = abs(ori_bd['x']['max'] - ori_bd['x']['min'])

    assert __is_match(nicv_ori.xs, nicv_nxt.xs, ref, total_count)
    assert __is_match(nicv_ori.ys, nicv_nxt.ys, ref, total_count)


def test_datatable_13C_DEPT135():
    nicv_ori, nicp_ori, jcamp_ori = __generated_jcamp_temp(
        __fixture_path(C13_DEPT135_dx)
    )
    nicv_nxt, nicp_nxt, jcamp_nxt = __generated_jcamp_temp(jcamp_ori.name)
    total_count = nicv_ori.ys.shape[0]
    ori_bd = nicv_ori.boundary
    ref = abs(ori_bd['x']['max'] - ori_bd['x']['min'])

    assert __is_match(nicv_ori.xs, nicv_nxt.xs, ref, total_count)
    assert __is_match(nicv_ori.ys, nicv_nxt.ys, ref, total_count)


def test_datatable_SVS_790A_13C_jdx():
    nicv_ori, nicp_ori, jcamp_ori = __generated_jcamp_temp(
        __fixture_path(SVS_790A_13C_jdx)
    )
    nicv_nxt, nicp_nxt, jcamp_nxt = __generated_jcamp_temp(jcamp_ori.name)
    total_count = nicv_ori.ys.shape[0]
    ori_bd = nicv_ori.boundary
    ref = abs(ori_bd['x']['max'] - ori_bd['x']['min'])

    assert __is_match(nicv_ori.xs, nicv_nxt.xs, ref, total_count)
    assert __is_match(nicv_ori.ys, nicv_nxt.ys, ref, total_count, 1/10000)

def test_datatable_JPK_948_jdx():
    nicv_ori, nicp_ori, jcamp_ori = __generated_jcamp_temp(
        __fixture_path(JPK_948_jdx)
    )
    nicv_nxt, nicp_nxt, jcamp_nxt = __generated_jcamp_temp(jcamp_ori.name)
    total_count = nicv_ori.ys.shape[0]
    ori_bd = nicv_ori.boundary
    ref = abs(ori_bd['x']['max'] - ori_bd['x']['min'])

    assert __is_match(nicv_ori.xs, nicv_nxt.xs, ref, total_count)
    assert __is_match(nicv_ori.ys, nicv_nxt.ys, ref, total_count)

def test_datatable_ms():
    mscv_ori, mscp_ori, jcamp_ori = __generated_jcamp_temp(
        __fixture_path(MS_dx)
    )
    mscv_nxt, mscp_nxt, jcamp_nxt = __generated_jcamp_temp(jcamp_ori.name)
    
    assert mscv_ori.datatables == mscv_nxt.datatables
    assert mscp_ori.core.datatables == mscp_nxt.core.datatables
