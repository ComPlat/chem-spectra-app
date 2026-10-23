import tempfile
import zipfile
from flask import current_app

from chem_spectra.lib.converter.share import UnconvertibleSpectrum

# Limits on what an uploaded archive may unpack to. Both are deployment
# settings (`MAX_ZIP_SIZE`, `MAX_ZIP_MEMBERS`, or `FLASK_`-prefixed in the
# environment); these apply when neither is set.
DEFAULT_MAX_ZIP_SIZE = 200      # MB uncompressed; BagIt LC/MS archives reach ~120
DEFAULT_MAX_ZIP_MEMBERS = 1000  # the largest Bruker or BagIt zip seen holds 38


def _limit(key, default):
    try:
        value = current_app.config.get(key)
    except RuntimeError:
        # outside an application context, e.g. a converter test
        value = None
    return default if value is None else value


def refuse_if_archive_too_large(archive):
    """Refuse an archive whose contents exceed the configured limits.

    Judged from the archive's directory, before anything is extracted.
    Python's zipfile stops reading a member at its declared size, so the
    declared sizes are the ones extraction will produce.
    """
    members = archive.infolist()
    max_members = _limit('MAX_ZIP_MEMBERS', DEFAULT_MAX_ZIP_MEMBERS)
    if len(members) > max_members:
        raise UnconvertibleSpectrum(
            'this archive holds {} files; the limit is {}.'.format(
                len(members), max_members))
    size_mb = sum(e.file_size for e in members) / (1024 * 1024)
    max_size_mb = _limit('MAX_ZIP_SIZE', DEFAULT_MAX_ZIP_SIZE)
    if size_mb > max_size_mb:
        raise UnconvertibleSpectrum(
            'this archive expands to {:.1f} MB; the limit is {} MB.'.format(
                size_mb, max_size_mb))


class FileContainer:
    def __init__(self, src=False):
        self.name = src and src.filename
        self.mimetype = src and src.mimetype
        self.bcore = src and src.stream.read()
        if src and zipfile.is_zipfile(src):
            with zipfile.ZipFile(src) as z:
                refuse_if_archive_too_large(z)

        self.core = self.bcore and self.bcore.decode('utf-8', errors='ignore')

    def from_str(self, core):
        self.core = core
        self.bcore = str.encode(core)
        return self

    def temp_file(self, prefix=None):
        suffix = '.{}'.format(self.name.split('.')[-1])
        tf = tempfile.NamedTemporaryFile(suffix=suffix)
        tf.write(self.bcore)
        tf.seek(0)
        return tf
