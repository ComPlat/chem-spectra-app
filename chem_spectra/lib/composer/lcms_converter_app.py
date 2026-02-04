import tempfile
from typing import List, Optional


class LCMSConverterAppComposer:
    def __init__(
        self,
        jcamp_files: List[tempfile.NamedTemporaryFile],
        image: Optional[tempfile.NamedTemporaryFile] = None,
    ):
        self.data = jcamp_files
        self._image = image

    def tf_img(self):
        return self._image

    def tf_jcamp(self):
        return self.data[0] if self.data else None
