import tempfile
import zipfile
from flask import current_app
import logging

DEFAULT_MAX_ZIP_SIZE = 100  # 100MB


class FileContainer:
    def __init__(self, src=False):
        self.name = src and src.filename
        self.mimetype = src and src.mimetype
        self.bcore = src and src.stream.read()
        if (self.mimetype in ['application/zip', 'application/octet-stream']) and zipfile.is_zipfile(src):
            with zipfile.ZipFile(src) as z:
                total_file_size = sum(e.file_size for e in z.infolist())
                try:
                    max_zip_size = current_app.config.get('MAX_ZIP_SIZE')
                    if max_zip_size is None:
                        max_zip_size = DEFAULT_MAX_ZIP_SIZE
                except Exception as exception:
                    max_zip_size = DEFAULT_MAX_ZIP_SIZE
                total_size_in_MB = total_file_size/(1024*1024)
                if total_size_in_MB > max_zip_size:
                    logger = logging.getLogger(__name__)
                    logger.setLevel(logging.ERROR)
                    logger.error(
                        f'This is maybe a zip bombs because its size is {total_size_in_MB} MB but maximum size of zip file is {max_zip_size} MB')
                    self.bcore = None

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
