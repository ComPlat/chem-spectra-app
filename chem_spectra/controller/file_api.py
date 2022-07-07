import base64
from flask import (
    Blueprint, request, jsonify, send_file, abort,
)

# from chem_spectra.controller.helper.settings import get_ip_white_list
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.controller.helper.share import (
    to_zip_response, extract_params, to_zip_bag_it_response
)
from chem_spectra.model.transformer import TransformerModel as TraModel
from chem_spectra.model.molecule import MoleculeModel
from chem_spectra.lib.converter.bagit.base import BagItBaseConverter


file_api = Blueprint('file_api', __name__)


@file_api.route('/api/v1/chemspectra/file/convert', methods=['POST'])
def chemspectra_file_convert():
    file = FileContainer(request.files['file'])
    molfile = FileContainer(request.files.get('molfile'))
    params = extract_params(request)
    if file:
        tf_jcamp, tf_img, tf_csv = TraModel(file, molfile=molfile, params=params).convert2jcamp_img()
        if not tf_jcamp:
            if isinstance(tf_img, BagItBaseConverter):
                list_jcamps = tf_img.get_base64_data()
                return jsonify(
                    status=True,
                    list_jcamps=list_jcamps,
                )
            else:
                abort(400)
        jcamp = base64.b64encode(tf_jcamp.read()).decode("utf-8")
        img = base64.b64encode(tf_img.read()).decode("utf-8")
        return jsonify(
            status=True,
            jcamp=jcamp,
            img=img
        )


@file_api.route('/api/v1/chemspectra/file/save', methods=['POST'])
def chemspectra_file_save():
    src = FileContainer(request.files['src'])
    request_files = request.files
    molfile = FileContainer(request.files.get('molfile'))
    filename = request.form.get('filename', default=None)
    params = extract_params(request)
    if 'dst_list' in request_files:
        request_dst = request_files.getlist('dst_list')
        dst_list = []
        for index in range(len(request_dst)):
            item = request_dst[index]
            dst = FileContainer(item)
            params['jcamp_idx'] = index
            if dst:  # and allowed_file(file):
                tm = TraModel(dst, molfile=molfile, params=params)
                tf_jcamp, tf_img, tf_csv = tm.convert2jcamp_img()
                if (tf_csv is not None and tf_csv != False):
                    tf_arr = [
                        src.temp_file(),
                        tf_jcamp,
                        tf_img,
                        tm.tf_predict(),
                        tf_csv,
                    ]
                else:
                    tf_arr = [
                        src.temp_file(),
                        tf_jcamp,
                        tf_img,
                        tm.tf_predict(),
                    ]
                dst_list.append(tf_arr)
        memory = to_zip_bag_it_response(dst_list, filename, src_idx=0)
        return send_file(
            memory,
            attachment_filename='spectrum.zip',
            as_attachment=True
        )
    else:
        dst = FileContainer(request.files['dst'])
        if dst:  # and allowed_file(file):
            tm = TraModel(dst, molfile=molfile, params=params)
            tf_jcamp, tf_img, tf_csv = tm.convert2jcamp_img()
            if (tf_csv is not None and tf_csv != False):
                tf_arr = [
                    src.temp_file(),
                    tf_jcamp,
                    tf_img,
                    tm.tf_predict(),
                    tf_csv,
                ]
            else:
                tf_arr = [
                    src.temp_file(),
                    tf_jcamp,
                    tf_img,
                    tm.tf_predict(),
                ]
            memory = to_zip_response(tf_arr, filename, src_idx=0)
            return send_file(
                memory,
                attachment_filename='spectrum.zip',
                as_attachment=True
            )


@file_api.route('/api/v1/chemspectra/file/refresh', methods=['POST'])
def chemspectra_file_refresh():
    # src = FileContainer(request.files['src'])
    request_files = request.files
    # dst = FileContainer(request.files['dst'])
    molfile = FileContainer(request.files.get('molfile'))
    # filename = request.form.get('filename', default=None)
    params = extract_params(request)
    if 'dst_list' in request_files:
        request_dst = request_files.getlist('dst_list')
        dst_list = []
        for index in range(len(request_dst)):
            item = request_dst[index]
            dst = FileContainer(item)
            params['jcamp_idx'] = index
            if dst:  # and allowed_file(file):
                tm = TraModel(dst, molfile=molfile, params=params)
                tf_jcamp, tf_img, tf_csv = tm.convert2jcamp_img()
                if (tf_csv is not None and tf_csv != False):
                    tf_arr = [
                        tf_jcamp,
                        tf_img,
                        tm.tf_predict(),
                        tf_csv
                    ]
                else:
                    tf_arr = [
                        tf_jcamp,
                        tf_img,
                        tm.tf_predict(),
                    ]
                dst_list.append(tf_arr)
        memory = to_zip_bag_it_response(dst_list, 'spectrum.zip')
        return send_file(
            memory,
            attachment_filename='spectrum.zip',
            as_attachment=True
        )
    else:
        dst = FileContainer(request.files['dst'])
        if dst:  # and allowed_file(file):
            tm = TraModel(dst, molfile=molfile, params=params)
            tf_jcamp, tf_img = tm.convert2jcamp_img()
            if not tf_jcamp:
                abort(400)
            jcamp = base64.b64encode(tf_jcamp.read()).decode("utf-8")
            img = base64.b64encode(tf_img.read()).decode("utf-8")
            return jsonify(
                status=True,
                jcamp=jcamp,
                img=img
            )


@file_api.route('/api/v1/chemspectra/molfile/convert', methods=['POST'])
def chemspectra_molfile_convert():
    molfile = FileContainer(request.files.get('molfile'))
    mm = MoleculeModel(molfile)
    return jsonify(
        status=True,
        smi=mm.smi,
        mass=mm.mass,
        svg=mm.svg,
    )
