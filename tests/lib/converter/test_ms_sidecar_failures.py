"""When msconvert does not deliver, say so — quickly, and with a reason.

Every failure of the msconvert sidecar used to arrive the same way: two
minutes of polling an empty directory, then

    TypeError: 'NoneType' object is not iterable   -- ms.py:296

because `__read_mz_ml` returned `(None, None, 0)` on expiry and
`__set_datatables` iterated it. The cause was discarded three times over — the
sidecar returns a status, the `/bin/docker` shim prints it and exits 0, and
`__run_cmd` ignored the exit code, stdout and stderr.

These tests pin the part that is fixable in this repository alone. Propagating
the sidecar's *HTTP* status needs a change to the shim, which lives in
chemotion-bakery; the shim already exits non-zero when it cannot reach the
service at all, and that is the failure these cover.
"""

import os
import subprocess as sbp

import pytest
from werkzeug.datastructures import FileStorage

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.lib.converter.ms import (
    MSConversionFailed,
    MSConverter,
    SIDECAR_TIMEOUT,
)

RAW = './tests/fixtures/source/ms/MS_ESI.RAW'
MZML = './tests/fixtures/source/ms/svs813f1.mzML'
PARAMS = {'mass': 230.079907196}


def _upload(path):
    handle = open(path, 'rb')
    return FileContainer(FileStorage(handle)), handle


def _fake_run(returncode=0, stdout='', stderr='', raises=None):
    def run(cmd, **kwargs):
        if raises is not None:
            raise raises
        return sbp.CompletedProcess(cmd, returncode, stdout, stderr)
    return run


# - - - the sidecar could not be reached - - -

def test_an_unreachable_sidecar_fails_immediately_and_says_so(monkeypatch):
    """The measured failure. The shim exits 1 because its `requests.post`
    raises outside its own try/except, so the exit code alone is enough to
    catch this without any change to the shim."""
    monkeypatch.setattr(
        sbp, 'run',
        _fake_run(returncode=1, stderr='ConnectionError: [Errno -2] Name or service not known'))

    upload, handle = _upload(RAW)
    try:
        with pytest.raises(MSConversionFailed) as caught:
            MSConverter(upload, PARAMS)
    finally:
        handle.close()

    message = str(caught.value)
    assert 'msconvert service' in message
    assert 'ConnectionError' in message, 'the sidecar\'s own words must survive'


def test_the_failure_is_not_the_old_nonetype_typeerror(monkeypatch):
    """What regressing would look like: `TypeError: 'NoneType' object is not
    iterable`, two minutes later, naming nothing."""
    monkeypatch.setattr(sbp, 'run', _fake_run(returncode=1, stderr='boom'))

    upload, handle = _upload(RAW)
    try:
        with pytest.raises(MSConversionFailed):
            MSConverter(upload, PARAMS)
    except TypeError:                                   # pragma: no cover
        pytest.fail('the NoneType failure is back')
    finally:
        handle.close()


def test_a_hanging_sidecar_is_bounded(monkeypatch):
    monkeypatch.setattr(
        sbp, 'run',
        _fake_run(raises=sbp.TimeoutExpired(cmd='docker', timeout=SIDECAR_TIMEOUT)))

    upload, handle = _upload(RAW)
    try:
        with pytest.raises(MSConversionFailed) as caught:
            MSConverter(upload, PARAMS)
    finally:
        handle.close()
    assert str(SIDECAR_TIMEOUT) in str(caught.value)


def test_our_timeout_leaves_room_for_the_sidecars_own():
    """`mscrunner.py` runs msconvert under `timeout=10` and answers with the
    result. Ours has to be longer, or we cut off the reply we asked for — when
    both were 10s they raced and the sidecar's account was usually lost."""
    assert SIDECAR_TIMEOUT > 10


# - - - the conversion ran but produced nothing - - -

def test_a_missing_mzml_is_reported_rather_than_iterated(monkeypatch):
    """Sidecar reports success, no file appears. Previously this polled for
    120s and then raised TypeError somewhere else entirely."""
    monkeypatch.setattr(sbp, 'run', _fake_run(returncode=0))
    monkeypatch.setattr('chem_spectra.lib.converter.ms.MZML_WAIT', 0.05)

    upload, handle = _upload(RAW)
    try:
        with pytest.raises(MSConversionFailed) as caught:
            MSConverter(upload, PARAMS)
    finally:
        handle.close()
    assert 'mzML' in str(caught.value)


# - - - cleanup - - -

def test_a_failed_conversion_leaves_nothing_behind(monkeypatch):
    """`__clean()` used to be the last statement of `__init__`, after the
    raise, so every failed upload leaked its hashed directory and the uploaded
    file onto the shared volume — the one the msconvert sidecar also mounts."""
    monkeypatch.setattr(sbp, 'run', _fake_run(returncode=1, stderr='nope'))

    captured = {}
    from chem_spectra.lib.converter import ms as ms_module
    original = ms_module.MSConverter._MSConverter__mk_dir

    def remember(self):
        target_dir, hash_str = original(self)
        captured['dir'] = target_dir
        return target_dir, hash_str

    monkeypatch.setattr(ms_module.MSConverter, '_MSConverter__mk_dir', remember)

    upload, handle = _upload(RAW)
    try:
        with pytest.raises(MSConversionFailed):
            MSConverter(upload, PARAMS)
    finally:
        handle.close()

    assert captured['dir'] is not None
    assert not os.path.exists(captured['dir'].absolute().as_posix()), \
        'the hashed directory survived a failed conversion'


# - - - the paths that never touch the sidecar - - -

def test_an_mzml_upload_is_unaffected():
    """mzML and mzXML never invoke msconvert, so none of the above applies to
    them. This is the guard that the changes did not reach that path."""
    upload, handle = _upload(MZML)
    try:
        converter = MSConverter(upload, PARAMS)
    finally:
        handle.close()
    assert converter.typ == 'MS'
    assert converter.spectra is not None
    assert len(converter.datatables) > 0
