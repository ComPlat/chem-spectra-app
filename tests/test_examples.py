import pytest
import io

from werkzeug.datastructures import FileStorage
from chem_spectra.lib.spectra.helper import (
    convert2jcamp_img, create_sp_carrier
)
from chem_spectra.lib.spectra.writer import (
    build_block_meta
)

from test_spectra_peaks import (
    __generated_peaks_meta, __target_peaks_meta,
)


params_0 = {
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# MNova
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def test_mnova_13C():
    target = 'mnova/MNOVA_SVS_13C.jdx'
    inherit_meta_content = __generated_peaks_meta(target, params_0)
    meta_target = __target_peaks_meta(target)
    assert inherit_meta_content == meta_target


def test_STM212_H():
    target = 'mnova/STM212_H.jcamp'
    inherit_meta_content = __generated_peaks_meta(target, params_0)
    meta_target = __target_peaks_meta(target)
    assert inherit_meta_content == meta_target


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# MNova pick-peaking
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def test_mnova_nonpp():
    target = 'mnova/mnova_nonpp.jcamp'
    inherit_meta_content = __generated_peaks_meta(target, params_0)
    meta_target = __target_peaks_meta(target)
    assert inherit_meta_content == meta_target


def test_mnova_haspp():
    target = 'mnova/mnova_haspp.jcamp'
    inherit_meta_content = __generated_peaks_meta(target, params_0)
    meta_target = __target_peaks_meta(target)
    assert inherit_meta_content == meta_target


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# Example
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def test_fdm_175():
    target = 'example/fdm_c175.dx'
    inherit_meta_content = __generated_peaks_meta(target, params_0)
    meta_target = __target_peaks_meta(target)
    assert inherit_meta_content == meta_target
