"""LCMS-specific Flask response builders.

Encapsulates the multi-file LCMS upload path served by
``/zip_jcamp_n_img`` so the generic transform_api stays focused on the
canonical (single-file/BagIt) flow.
"""

import json
from typing import Iterable, List, Optional

from flask import Response, make_response, send_file

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.controller.helper.share import to_zip_bag_it_response
from chem_spectra.lib.composer.lcms_converter_app import LCMSConverterAppComposer


def _materialize_uploads(uploads: Iterable) -> List:
    """Convert raw Flask file uploads into closeable temporary files."""
    temp_files = []
    for item in uploads:
        container = FileContainer(item)
        if container and container.bcore:
            temp_files.append(container.temp_file())
    return temp_files


def _close_silently(handles: Iterable) -> None:
    for handle in handles:
        try:
            handle.close()
        except Exception:
            pass


def build_lcms_zip_response_from_uploads(
    uploads: Iterable,
    params: Optional[dict],
    download_name: str = 'spectrum.zip',
) -> Response:
    """Bundle a batch of LCMS JCAMP uploads into a BagIt-style zip response.

    Mirrors what ``BagItBaseConverter`` produces for a BagIt LCMS dataset, but
    starts from individual files instead of an archive.
    """
    temp_files = _materialize_uploads(uploads)
    try:
        composer = LCMSConverterAppComposer(temp_files, None, params)

        updated_jcamp = composer.tf_jcamp()
        if updated_jcamp and updated_jcamp not in temp_files:
            temp_files.append(updated_jcamp)

        tf_img = composer.tf_img()
        list_jcamps = composer.data or temp_files

        dst_list = []
        for idx, tf_jcamp_file in enumerate(list_jcamps):
            if idx == 0 and tf_img is not None:
                dst_list.append([tf_jcamp_file, tf_img])
            else:
                dst_list.append([tf_jcamp_file])

        memory = to_zip_bag_it_response(dst_list)
        rsp = make_response(send_file(
            memory,
            download_name=download_name,
            as_attachment=True,
        ))
        rsp.headers['X-Extra-Info-JSON'] = json.dumps(
            {'spc_type': 'hplc', 'invalid_molfile': False}
        )
        return rsp
    finally:
        _close_silently(temp_files)
