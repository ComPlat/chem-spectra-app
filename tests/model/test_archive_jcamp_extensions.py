"""A JCAMP file inside an archive is read whichever spelling it has.

`.jdx`, `.dx` and `.jcamp` are one format. A `.dx` or `.jcamp` file
converted when uploaded on its own, but the archive readers looked for
`.jdx` alone: in a BagIt the member was skipped without a word, and a plain
zip holding only those spellings was not recognised as JCAMP at all.
"""

import io
import zipfile

import pytest
from werkzeug.datastructures import FileStorage

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.lib.converter.bagit.base import BagItBaseConverter
from chem_spectra.lib.converter.share import is_jcamp_filename
from chem_spectra.model.transformer import TransformerModel

CV_BAGIT = './tests/fixtures/source/bagit/cv/File053_BagIt.zip'
IR = './tests/fixtures/source/IR.dx'
UVVIS = './tests/fixtures/source/JPK-948.jdx'


@pytest.mark.parametrize('name, expected', [
    ('a.jdx', True), ('a.dx', True), ('a.jcamp', True),
    ('A.JDX', True), ('a.peak.jdx', True),
    ('a.txt', False), ('a.json', False), ('dx', False), ('a.jdx.bak', False),
])
def test_is_jcamp_filename(name, expected):
    assert is_jcamp_filename(name) is expected


def test_a_bagit_reads_every_jcamp_spelling(tmp_path):
    with zipfile.ZipFile(CV_BAGIT) as archive:
        archive.extractall(tmp_path)
    data = tmp_path / 'data'
    (data / 'table_02.jdx').rename(data / 'table_02.dx')
    (data / 'table_03.jdx').rename(data / 'table_03.jcamp')

    converter = BagItBaseConverter(str(tmp_path))

    assert converter.archive_entry_stems == ['table_01', 'table_02', 'table_03']
    assert len(converter.data) == 3
    assert len(converter.images) == 3


def _zip(members):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, 'w') as archive:
        for name, path in members.items():
            with open(path, 'rb') as handle:
                archive.writestr(name, handle.read())
    buffer.seek(0)
    return FileContainer(FileStorage(stream=buffer, filename='spectra.zip',
                                     content_type='application/zip'))


def test_a_plain_zip_without_any_jdx_is_recognised():
    file = _zip({'ir.dx': IR, 'uvvis.jcamp': UVVIS, 'notes.txt': IR})
    converter, composer, _ = TransformerModel(
        file, params={'fname': 'spectra.zip', 'ext': 'zip'}).zip2cvp()

    assert isinstance(converter, BagItBaseConverter)
    assert composer is converter
    # notes.txt is not a JCAMP name, whatever it holds
    assert converter.archive_entry_stems == ['ir', 'uvvis']
