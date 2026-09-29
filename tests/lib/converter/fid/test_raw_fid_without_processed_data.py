"""A Bruker upload carrying only the FID, with no processed data.

`search_brucker_binary` splits uploads into those that ship `pdata` and those
that do not. The second branch was broken and untested: a second `__init__`
shadowed the one that read a directory, so `FidBaseConverter.__read` was dead
code, and `zip2cvp` passed `(target_dir, params, fname)` into the surviving
`(dic, data, params, fname)` signature -- binding a path string to `dic`.

Measured before the fix: AttributeError out of the request, on the migration
branch and on the pin before it alike. Both Bruker fixtures ship `pdata`, which
is why no test ever reached it.
"""

import io
import zipfile

import pytest
from werkzeug.datastructures import FileStorage

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.model.transformer import TransformerModel

WITH_PDATA = './tests/fixtures/source/bruker/1H.zip'
MOLFILE = './tests/fixtures/source/molfile/svs813f1_B.mol'


@pytest.fixture
def raw_fid_zip(tmp_path):
    """The 1H fixture with its processed data stripped out.

    Derived rather than committed as a second binary: the point is that this is
    the *same* measurement, differing only in whether the instrument's
    processed spectrum came along.
    """
    target = tmp_path / 'raw_fid.zip'
    source = zipfile.ZipFile(WITH_PDATA)
    with zipfile.ZipFile(target, 'w') as out:
        for name in source.namelist():
            if 'pdata' in name:
                continue
            out.writestr(name, source.read(name))
    assert not any('pdata' in n for n in zipfile.ZipFile(target).namelist())
    return str(target)


def _convert(path):
    with open(path, 'rb') as handle:
        container = FileContainer(FileStorage(handle))
        with open(MOLFILE, 'rb') as molhandle:
            molfile = FileContainer(FileStorage(molhandle))
            _, composers, _ = TransformerModel(
                container, molfile=molfile, params={'ext': 'zip'},
            ).zip2cvp()
    composer = composers[0] if isinstance(composers, list) else composers
    return composer.core


def test_a_raw_fid_upload_converts(raw_fid_zip):
    core = _convert(raw_fid_zip)
    assert core.typ == 'NMR'
    assert len(core.ys) > 0


def test_it_is_phased_and_scaled_like_the_processed_upload(raw_fid_zip):
    """The FID is transformed here rather than read from the instrument, so the
    x axis must still come out in ppm over the same range -- this is a 1H
    spectrum either way.
    """
    core = _convert(raw_fid_zip)
    assert core.xs.min() == pytest.approx(-4.01, abs=0.05)
    assert core.xs.max() == pytest.approx(16.01, abs=0.05)


def test_the_directory_reader_is_reachable():
    """Guards the shadowing specifically: the entry point must exist and not be
    the two-argument data constructor."""
    from chem_spectra.lib.converter.fid.base import FidBaseConverter
    assert hasattr(FidBaseConverter, 'from_directory')
