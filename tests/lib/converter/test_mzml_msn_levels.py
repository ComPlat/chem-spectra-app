"""An mzML with MS^4 or deeper spectra was refused as unreadable.

pymzml looks up a measured precision per MS level and knows levels 0-3 by
default, so iterating over an MS^4 spectrum raised KeyError: 4 inside the
reader, and the upload came back as "could not be read as mzML: KeyError: 4".
The test file is the repo's own svs813f1.mzML with some spectra relabelled
to deeper MS levels; the edit keeps every byte offset.
"""
import io

from werkzeug.datastructures import FileStorage

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.lib.composer.ms import MSComposer
from chem_spectra.lib.converter.ms import MSConverter

SOURCE = './tests/fixtures/source/ms/svs813f1.mzML'
LEVEL_1 = b'name="ms level" value="1"'


def _with_levels(levels):
    with open(SOURCE, 'rb') as f:
        content = f.read()
    for level in levels:
        content = content.replace(
            LEVEL_1, b'name="ms level" value="%d"' % level, 1)
    return content


def test_spectra_beyond_ms3_are_read():
    content = _with_levels([2, 3, 4, 5])
    upload = FileStorage(io.BytesIO(content), filename='msn.mzML')

    mscv = MSConverter(FileContainer(upload), {'mass': 230.079907196})
    lines = MSComposer(mscv).tf_jcamp().read().decode('utf-8').split('\n')

    assert '##$CSSCANCOUNT=24' in lines
