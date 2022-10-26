from cmath import log
from crypt import methods
import json
import collections.abc
from multiprocessing.dummy import Array
from flask import (
    Blueprint, request, send_file, make_response,
)
# from chem_spectra.controller.helper.settings import get_ip_white_list
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.controller.helper.share import (
    to_zip_response, extract_params, to_zip_bag_it_response
)

from chem_spectra.model.transformer import TransformerModel as TraModel
from chem_spectra.lib.converter.bagit.base import BagItBaseConverter


trans_api = Blueprint('transform_api', __name__)


@trans_api.before_app_request
def filter_remote_ip():
    pass
    # trusted_servers = get_ip_white_list()
    # if request.remote_addr not in trusted_servers:
    #     abort(403)


@trans_api.route('/zip_jcamp_n_img', methods=['POST'])
def zip_jcamp_n_img():
    file = FileContainer(request.files['file'])
    molfile = FileContainer(request.files.get('molfile'))
    params = extract_params(request)
    if file:  # and allowed_file(file):
        cmpsr = TraModel(file, molfile=molfile, params=params).to_composer()
        if ((type(cmpsr) is dict) and "invalid_molfile" in cmpsr):
            return json.dumps(cmpsr)

        if isinstance(cmpsr, BagItBaseConverter):
            # check if composered model is in BagIt format
            list_jcamps, list_images, list_csv = cmpsr.data, cmpsr.images, cmpsr.list_csv
            dst_list = []
            for idx in range(len(list_jcamps)):
                tf_jcamp = list_jcamps[idx]
                tf_img = list_images[idx]
                tf_csv = list_csv[idx]
                tf_arr = [tf_jcamp, tf_img, tf_csv]
                dst_list.append(tf_arr)

            memory = to_zip_bag_it_response(dst_list)
            rsp = make_response(
                send_file(
                    memory,
                    attachment_filename='spectrum.zip',
                    as_attachment=True
                )
            )
            rsp.headers['X-Extra-Info-JSON'] = json.dumps({'spc_type': 'bagit'})
        elif isinstance(cmpsr, collections.abc.Sequence):
            dst_list = []
            spc_type = ''
            for composer in cmpsr:
                spc_type = composer.core.ncl if composer.core.typ == 'NMR' else composer.core.typ
                tf_jcamp, tf_img = composer.tf_jcamp(), composer.tf_img()
                tf_arr = [tf_jcamp, tf_img]
                dst_list.extend(tf_arr)
                
            memory = to_zip_response(dst_list)
            rsp = make_response(
                send_file(
                    memory,
                    attachment_filename='spectrum.zip',
                    as_attachment=True
                )
            )
            rsp.headers['X-Extra-Info-JSON'] = json.dumps({'spc_type': spc_type})
        else:
            tf_jcamp, tf_img, tf_csv = cmpsr.tf_jcamp(), cmpsr.tf_img(), cmpsr.tf_csv()
            spc_type = cmpsr.core.ncl if cmpsr.core.typ == 'NMR' else cmpsr.core.typ
            if (tf_csv is not None and tf_csv != False):
                memory = to_zip_response([tf_jcamp, tf_img, tf_csv])
            else:
                memory = to_zip_response([tf_jcamp, tf_img])
            rsp = make_response(
                send_file(
                    memory,
                    attachment_filename='spectrum.zip',
                    as_attachment=True
                )
            )
            rsp.headers['X-Extra-Info-JSON'] = json.dumps({'spc_type': spc_type})
        return rsp


@trans_api.route('/zip_jcamp', methods=['POST'])
def zip_jcamp():
    file = FileContainer(request.files['file'])
    molfile = FileContainer(request.files.get('molfile'))
    params = extract_params(request)
    if file:  # and allowed_file(file):
        tf_jcamp = TraModel(file, molfile=molfile, params=params).convert2jcamp()
        if isinstance(tf_jcamp, BagItBaseConverter):
            # check if composered model is in BagIt format
            list_jcamps = tf_jcamp.data
            dst_list = []
            for jcamp in list_jcamps:
                tf_arr = [jcamp]
                dst_list.append(tf_arr)

            memory = to_zip_bag_it_response(dst_list)
        else:
            memory = to_zip_response([tf_jcamp])
        return send_file(
            memory,
            attachment_filename='spectrum.zip',
            as_attachment=True
        )


@trans_api.route('/zip_image', methods=['POST'])
def zip_image():
    file = FileContainer(request.files['file'])
    molfile = FileContainer(request.files.get('molfile'))
    params = extract_params(request)
    if file:  # and allowed_file(file):
        tf_img = TraModel(file, molfile=molfile, params=params).convert2img()
        if isinstance(tf_img, BagItBaseConverter):
            # check if composered model is in BagIt format
            list_images = tf_img.images
            dst_list = []
            for img in list_images:
                tf_arr = [img]
                dst_list.append(tf_arr)

            memory = to_zip_bag_it_response(dst_list)
        else:
            memory = to_zip_response([tf_img])
        return send_file(
            memory,
            attachment_filename='spectrum.zip',
            as_attachment=True
        )


@trans_api.route('/jcamp', methods=['POST'])
def jcamp():
    file = FileContainer(request.files['file'])
    molfile = FileContainer(request.files.get('molfile'))
    params = extract_params(request)
    if file:  # and allowed_file(file):
        tf_jcamp = TraModel(file, molfile=molfile, params=params).convert2jcamp()
        return send_file(
            tf_jcamp,
            attachment_filename='spectrum.jdx',
            as_attachment=True
        )


@trans_api.route('/image', methods=['POST'])
def image():
    file = FileContainer(request.files['file'])
    molfile = FileContainer(request.files.get('molfile'))
    params = extract_params(request)
    if file:  # and allowed_file(file):
        tf_img = TraModel(file, molfile=molfile, params=params).convert2img()
        return send_file(
            tf_img,
            attachment_filename='spectrum.png',
            as_attachment=True,
            mimetype='image/png'
        )

@trans_api.route('/nmrium', methods=['POST'])
def nmrium():
    nmriumFile = FileContainer(request.files['file'])
    if nmriumFile:
        # print(nmriumFile.bcore)
        transformModel = TraModel(file=nmriumFile)
        transformedData = transformModel.tf_nmrium()
        return send_file(
            transformedData,
            attachment_filename='spectrum.jdx',
            as_attachment=True
        )
