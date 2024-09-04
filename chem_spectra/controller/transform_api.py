from cmath import log
from crypt import methods
import json
import collections.abc
from multiprocessing.dummy import Array
from flask import (
    Blueprint, request, send_file, make_response, abort
)
# from chem_spectra.controller.helper.settings import get_ip_white_list
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.controller.helper.share import (
    to_zip_response, extract_params, to_zip_bag_it_response
)

from chem_spectra.model.transformer import TransformerModel as TraModel
from chem_spectra.lib.composer.lcms import LCMSComposer
from chem_spectra.lib.converter.bagit.base import BagItBaseConverter
from chem_spectra.model.molecule import MoleculeModel


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
        transform_model = TraModel(file, molfile=molfile, params=params)
        cmpsr, invalid_molfile = transform_model.to_composer()
        if (not cmpsr):
            abort(403)

        if isinstance(cmpsr, BagItBaseConverter):
            # check if composered model is in BagIt format
            list_jcamps, list_images, list_csv, combined_image = cmpsr.data, cmpsr.images, cmpsr.list_csv, cmpsr.combined_image
            dst_list = []
            for idx in range(len(list_jcamps)):
                tf_jcamp = list_jcamps[idx]
                tf_img = list_images[idx]
                tf_csv = list_csv[idx]
                tf_arr = [tf_jcamp, tf_img, tf_csv]
                dst_list.append(tf_arr)
                
            if combined_image is not None:
                dst_list.append(combined_image)

            memory = to_zip_bag_it_response(dst_list)
            rsp = make_response(
                send_file(
                    memory,
                    download_name='spectrum.zip',
                    as_attachment=True
                )
            )
            rsp.headers['X-Extra-Info-JSON'] = json.dumps({'spc_type': 'bagit', 'invalid_molfile': invalid_molfile})
        elif isinstance(cmpsr, LCMSComposer):
            # check if composered model is hplc ms
            list_jcamps = cmpsr.data
            dst_list = []
            for idx in range(len(list_jcamps)):
                tf_jcamp = list_jcamps[idx]
                tf_arr = [tf_jcamp]
                dst_list.append(tf_arr)

            memory = to_zip_bag_it_response(dst_list)
            rsp = make_response(
                send_file(
                    memory,
                    download_name='spectrum.zip',
                    as_attachment=True
                )
            )
            rsp.headers['X-Extra-Info-JSON'] = json.dumps({'spc_type': 'hplc', 'invalid_molfile': invalid_molfile})
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
                    download_name='spectrum.zip',
                    as_attachment=True
                )
            )
            rsp.headers['X-Extra-Info-JSON'] = json.dumps({'spc_type': spc_type, 'invalid_molfile': invalid_molfile})
        else:
            tf_jcamp, tf_img, tf_csv = cmpsr.tf_jcamp(), cmpsr.tf_img(), cmpsr.tf_csv()
            tf_nmrium = None
            try:
                molecule_model = MoleculeModel(molfile, cmpsr.core.ncl, decorate=False)
                tf_nmrium = cmpsr.generate_nmrium(molfile_data=molecule_model.moltxt)
            except Exception:
                pass

            spc_type = cmpsr.core.ncl if cmpsr.core.typ == 'NMR' else cmpsr.core.typ
            if (tf_csv is not None and tf_csv != False):
                memory = to_zip_response([tf_jcamp, tf_img, tf_csv])
            elif (tf_nmrium is not None):
                memory = to_zip_response([tf_jcamp, tf_img, tf_nmrium])
            else:
                memory = to_zip_response([tf_jcamp, tf_img])
            rsp = make_response(
                send_file(
                    memory,
                    download_name='spectrum.zip',
                    as_attachment=True
                )
            )
            rsp.headers['X-Extra-Info-JSON'] = json.dumps({'spc_type': spc_type, 'invalid_molfile': invalid_molfile})

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
            download_name='spectrum.zip',
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
            download_name='spectrum.zip',
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
            download_name='spectrum.jdx',
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
            download_name='spectrum.png',
            as_attachment=True,
            mimetype='image/png'
        )

@trans_api.route('/nmrium', methods=['POST'])
def nmrium():
    nmriumFile = FileContainer(request.files['file'])
    if nmriumFile:
        transformModel = TraModel(file=nmriumFile)
        transformedData = transformModel.tf_nmrium()
        if transformedData is None:
            abort(404)
        return send_file(
            transformedData,
            download_name='spectrum.jdx',
            as_attachment=True
        )

@trans_api.route('/combine_images', methods=['POST'])
def combine_images():
    request_files = request.files
    if (not request_files):
        abort(400)
        
    list_files = []
    for file in request_files.getlist('files[]'):
        file_container = FileContainer(file)
        list_files.append(file_container)
    
    params = extract_params(request)
    extras = request.form.get('extras', default=None)

    transform_model = TraModel(None, params=params, multiple_files=list_files)
    tf_combine = transform_model.tf_combine(list_file_names=params['list_file_names'], extraParams=extras)
    if (not tf_combine):
        abort(400)
    
    memory = to_zip_response([tf_combine])
    return send_file(
        memory,
        download_name='spectrum.zip',
        as_attachment=True
    )
