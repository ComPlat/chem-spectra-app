from test_spectra_peaks import (
    __generated_peaks_meta, __target_peaks_meta,
)

target_dir = './tests/fixtures/'
source_dir = 'source/'
result_dir = 'result/'

IR_dx = 'IR.dx'
H1_dx = '1H.dx'
C13_CPD_dx = '13C-CPD.dx'
C13_DEPT135_dx = '13C-DEPT135.dx'

meta_IR_dx = 'meta_IR'
meta_H1_dx = 'meta_1H'
meta_C13_CPD_dx = 'meta_13C-CPD'
meta_C13_DEPT135_dx = 'meta_13C-DEPT135'


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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# given empty params
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
empty_params = {
    'integration': '{"stack":[],"refArea":1,"refFactor":1,"shift":0}',
    'multiplicity': '{"stack":[],"shift":0,"smExtext":false}',
}

def test_empty_params_meta_IR_dx():
    empty_meta_content = __generated_peaks_meta(IR_dx, empty_params)
    meta_target = __target_peaks_meta('im/im_' + meta_IR_dx)
    assert empty_meta_content == meta_target

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# given Integration/Multiplicity params
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
im_params = {
    'integration': '{"stack":[{"xL":0.5,"xU":1.3,"area":2.123},{"xL":2.6902344414431028,"xU":3.3735489740043914,"area":1.0728174301051716},{"xL":7.094987150286099,"xU":8.125215512171152,"area":0.25925501835212633}],"refArea":1.0728174301051716,"refFactor":"2","shift":0}',
    'multiplicity': '{"stack":[{"peaks":[{"x":3.2293749183216156,"y":811992.63},{"x":2.9536338778941817,"y":1338.16},{"x":2.701360160056788,"y":563.02}],"xExtent":{"xL":2.6902344414431028,"xU":3.3735489740043914},"yExtent":{"yL":-346913.762766988,"yU":1323092.2668168638},"mpyType":"t","js":[]},{"peaks":[{"x":7.87517143721676,"y":218025.88}],"xExtent":{"xL":7.094987150286099,"xU":8.125215512171152},"yExtent":{"yL":-396616.06022033,"yU":865828.6929629277},"mpyType":"m","js":[1.1,2.2,3.3]}],"shift":0,"smExtext":{"xL":2.6902344414431028,"xU":3.3735489740043914}}',
}

def test_im_params_meta_1H():
    im_meta_content = __generated_peaks_meta(H1_dx, im_params)
    meta_target = __target_peaks_meta('im/im_' + meta_H1_dx)
    assert im_meta_content == meta_target


def test_im_params_meta_13C_DEPT135():
    im_meta_content = __generated_peaks_meta(C13_DEPT135_dx, im_params)
    meta_target_A = __target_peaks_meta('im/im_' + meta_C13_DEPT135_dx)
    assert im_meta_content[:1500] == meta_target_A[:1500]
