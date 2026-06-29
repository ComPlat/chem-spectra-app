import json
from typing import Iterable, List, Optional

from flask import Response, make_response, send_file

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.controller.helper.share import to_zip_flat_bagit_response
from chem_spectra.lib.composer.lcms_converter_app import LCMSConverterAppComposer
from chem_spectra.lib.converter.bagit.lcms_builder import classify_lcms_stems


def _materialize_uploads(uploads: Iterable) -> List:
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
    upload_list = list(uploads)
    temp_files = _materialize_uploads(upload_list)
    try:
        composer = LCMSConverterAppComposer(temp_files, None, params)
        composer.tf_jcamp()
        preview = composer.tf_img()

        list_jcamps = composer.data or temp_files
        entry_stems = classify_lcms_stems([f.name for f in list_jcamps])
        if len(entry_stems) != len(list_jcamps):
            entry_stems = [str(i) for i in range(len(list_jcamps))]

        preview_idx = next(
            (i for i, s in enumerate(entry_stems) if s == 'lcms_uvvis'),
            0,
        )
        dst_list = [[f] for f in list_jcamps]
        if preview is not None and dst_list:
            dst_list[preview_idx].append(preview)

        archive_fname = next(
            (getattr(u, 'filename', None) for u in upload_list if getattr(u, 'filename', None)),
            None,
        ) or 'spectrum'

        memory = to_zip_flat_bagit_response(dst_list, archive_fname, entry_stems)
        rsp = make_response(send_file(
            memory,
            download_name=download_name,
            as_attachment=True,
        ))
        rsp.headers['X-Extra-Info-JSON'] = json.dumps(
            {'spc_type': 'lcms', 'invalid_molfile': False}
        )
        return rsp
    finally:
        _close_silently(temp_files)
