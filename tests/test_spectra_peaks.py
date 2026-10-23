from werkzeug.datastructures import FileStorage
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.model.transformer import TransformerModel as TraModel
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter


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

separator_e = '$$ === CHEMSPECTRA PEAK TABLE EDIT ==='
separator_a = '$$ === CHEMSPECTRA PEAK TABLE AUTO ==='
separator = '$$ === CHEMSPECTRA INTEGRALS AND MULTIPLETS ==='

def __generated_peaks_meta(orig_filename, params=False):
    with open(target_dir + source_dir + orig_filename, 'rb') as f:
        file = FileContainer(FileStorage(f))
        molfile = FileContainer(FileStorage(None))
        _, tcp, _ = TraModel(file, molfile, params).jcamp2cvp()

    parts = ''.join(tcp.meta).split(separator)
    assert len(parts) == 2

    meta_content = parts[1]
    return meta_content


def __generated_core(orig_filename, params=False):
    with open(target_dir + source_dir + orig_filename, 'rb') as f:
        file = FileContainer(FileStorage(f))
        molfile = FileContainer(FileStorage(None))
        nicv, _, _ = TraModel(file, molfile, params).jcamp2cvp()
    return nicv


def __target_peaks_meta(filename):
    meta_target = open(target_dir + result_dir + filename, 'rb')
    return meta_target.read().decode('utf-8', errors='ignore')


def __generated_composer(orig_filename, params=False):
    with open(target_dir + source_dir + orig_filename, 'rb') as f:
        file = FileContainer(FileStorage(f))
        molfile = FileContainer(FileStorage(None))
        _, composer, _ = TraModel(file, molfile, params).jcamp2cvp()
    return composer


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


def test_auto_peaks_exclude_1h_solvent_peaks():
    core = __generated_core(H1_dx)
    assert core.auto_peaks is not None
    solvent_auto_peaks = [peak_x for peak_x in core.auto_peaks['x'] if 2.48 < peak_x < 2.52]
    solvent_edit_peaks = [peak_x for peak_x in core.edit_peaks['x'] if 2.48 < peak_x < 2.52]
    assert len(solvent_auto_peaks) == 1
    assert len(solvent_edit_peaks) == 1
    assert abs(solvent_auto_peaks[0] - 2.50) < 0.01
    assert abs(solvent_edit_peaks[0] - 2.50) < 0.01


def test_meta_13C_CPD_dx():
    meta_content = __generated_peaks_meta(C13_CPD_dx)
    meta_target = __target_peaks_meta(meta_C13_CPD_dx)
    assert meta_content == meta_target


def test_meta_13C_DEPT135():
    meta_content = __generated_peaks_meta(C13_DEPT135_dx)
    meta_target = __target_peaks_meta(meta_C13_DEPT135_dx)
    assert meta_content[:1500] == meta_target[:1500]


def test_meta_SVS_790A_13C_jdx():
    meta_content = __generated_peaks_meta(SVS_790A_13C_jdx)
    meta_target = __target_peaks_meta(meta_SVS_790A_13C_jdx)
    assert meta_content == meta_target


def test_auto_peaks_exclude_13c_solvent_peaks():
    core = __generated_core(SVS_790A_13C_jdx)
    assert core.auto_peaks is not None
    # the chloroform triplet, +-32 Hz at 100 MHz; the window (74.16-80.16)
    # also holds sample peaks at 79.3 and 79.8, which stay
    solvent_auto_peaks = [peak_x for peak_x in core.auto_peaks['x'] if 76.4 < peak_x < 77.9]
    solvent_edit_peaks = [peak_x for peak_x in core.edit_peaks['x'] if 76.4 < peak_x < 77.9]
    assert len(solvent_auto_peaks) == 1
    assert len(solvent_edit_peaks) == 1
    assert abs(solvent_auto_peaks[0] - 77.16) < 0.2
    assert abs(solvent_edit_peaks[0] - 77.16) < 0.2


def test_solvent_filter_keeps_most_intense_peak_in_range():
    converter = JcampTechniqueConverter.__new__(JcampTechniqueConverter)
    converter.solv_peaks = [(74.16, 80.16)]
    converter.ncl = '13C'

    peaks = [
        {'x': 77.02, 'y': 90.0},
        {'x': 77.24, 'y': 100.0},
        {'x': 77.12, 'y': 80.0},
        {'x': 120.0, 'y': 70.0},
    ]

    filtered_peaks = converter._JcampTechniqueConverter__filter_solvent_peaks(peaks)
    kept_solvent_peaks = [peak for peak in filtered_peaks if 74.16 < peak['x'] < 80.16]

    assert len(kept_solvent_peaks) == 1
    assert kept_solvent_peaks[0]['x'] == 77.24
    assert kept_solvent_peaks[0]['y'] == 100.0



def test_solvent_filter_keeps_a_sample_peak_elsewhere_in_the_window():
    # the 13C chloroform window is +-3 ppm; only the triplet (+-32 Hz at
    # 100 MHz) is the solvent, not a sample signal 2.3 ppm away
    converter = JcampTechniqueConverter.__new__(JcampTechniqueConverter)
    converter.solv_peaks = [(74.16, 80.16)]
    converter.ncl = '13C'
    converter.obs_freq = 100.6

    peaks = [
        {'x': 77.00, 'y': 100.0},
        {'x': 77.32, 'y': 95.0},
        {'x': 76.68, 'y': 94.0},
        {'x': 79.30, 'y': 40.0},
    ]

    filtered = converter._JcampTechniqueConverter__filter_solvent_peaks(peaks)

    assert sorted(peak['x'] for peak in filtered) == [77.00, 79.30]

def test_meta_MS_jdx():
    pass
    # meta_content = __generated_peaks_meta('MS.dx')
    # import pdb; pdb.set_trace()
    # meta_target = __target_peaks_meta(meta_SVS_790A_13C_jdx)
    # assert meta_content == meta_target


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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# explicitly cleared peaks (peaks_str='') must not fall back
# to auto-detection nor keep any previously saved edit peaks
# https://github.com/ComPlat/chem-spectra-app/issues/288
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
clear_peaks_params = {
    'peaks_str': '',
}


def test_clear_peaks_removes_previously_saved_edit_peaks():
    meta_content = __generated_peaks_meta(
        'edit/edit_' + H1_dx, clear_peaks_params
    )
    edit_block = meta_content.split(separator_a)[0]
    assert '##NPOINTS=0\n' in edit_block
    assert '##PEAKTABLE= (XY..XY)\n##END=' in edit_block


def test_clear_peaks_does_not_auto_detect():
    meta_content = __generated_peaks_meta(H1_dx, clear_peaks_params)
    edit_block = meta_content.split(separator_a)[0]
    assert '##NPOINTS=0\n' in edit_block
    assert '##PEAKTABLE= (XY..XY)\n##END=' in edit_block

"""
def test_params_meta_1H():
    ps_meta_content = __generated_peaks_meta(H1_dx, params_1)
    meta_target = __target_peaks_meta('ps/ps_' + meta_H1_dx)
    assert ps_meta_content == meta_target


def test_params_meta_13C_DEPT135():
    ps_meta_content = __generated_peaks_meta(C13_DEPT135_dx, params_1)
    meta_target = __target_peaks_meta('ps/ps_' + meta_C13_DEPT135_dx)
    assert ps_meta_content[:1500] == meta_target[:1500]


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
    auto_meta_content = __generated_peaks_meta(
        'auto/auto_' + C13_DEPT135_dx, params_1
    )
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


def test_edit_params_meta_1H():
    auto_meta_content = __generated_peaks_meta('edit/edit_' + H1_dx, params_1)
    meta_target = __target_peaks_meta('auto/auto_ps_' + meta_H1_dx)
    assert auto_meta_content == meta_target


def test_edit_params_auto_meta_13C_DEPT135():
    auto_meta_content = __generated_peaks_meta(
        'edit/edit_' + C13_DEPT135_dx, params_1
        )
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
    assert empty_meta_content[:1500] == meta_target[:1500]


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
    inherit_meta_content = __generated_peaks_meta(
        'edit/inherit_' + H1_dx, params_3
    )
    meta_target = __target_peaks_meta('edit/inherit_' + meta_H1_dx)
    assert inherit_meta_content == meta_target
"""


def _acetone_converter(params, dic):
    # acetone-d6 13C: the CD3 septet around 29.84 and the carbonyl at 206.26,
    # plus one sample peak outside both windows
    converter = JcampTechniqueConverter.__new__(JcampTechniqueConverter)
    converter.ncl = '13C'
    converter.obs_freq = 100.6
    converter.solv_peaks = [(27.0, 33.0), (203.7, 209.7)]
    converter.params = params
    converter.dic = dic
    step = 19.5 / 100.6
    septet = [{'x': 29.84 + k * step, 'y': w}
              for k, w in zip(range(-3, 4), [1.0, 3.0, 6.0, 7.0, 6.0, 3.0, 1.0])]
    carbonyl = [{'x': 206.26, 'y': 5.0}]
    sample = [{'x': 128.4, 'y': 9.0}]
    return converter, septet + carbonyl + sample


def _filtered(converter, peaks):
    kept = converter._JcampTechniqueConverter__filter_solvent_peaks(peaks)
    return sorted(round(peak['x'], 2) for peak in kept)


def test_only_the_reference_solvent_signal_keeps_a_line():
    # the reference is the CD3 septet: its line stays for the reference, and
    # the carbonyl line, which nothing else would remove, goes
    converter, peaks = _acetone_converter(
        {'ref_name': None}, {'$CSSOLVENTVALUE': ['29.84']})
    assert _filtered(converter, peaks) == [29.84, 128.4]


def test_the_requested_reference_signal_is_the_one_kept():
    converter, peaks = _acetone_converter(
        {'ref_name': 'Acetone-d6 (broad)', 'ref_value': '206.26'},
        {'$CSSOLVENTVALUE': ['29.84']})
    assert _filtered(converter, peaks) == [128.4, 206.26]


def test_without_a_reference_value_the_first_window_is_the_reference():
    converter, peaks = _acetone_converter({'ref_name': None}, {})
    assert _filtered(converter, peaks) == [29.84, 128.4]
