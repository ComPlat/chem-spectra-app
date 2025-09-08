import nmrglue as ng
import numpy as np  # noqa: F401
import os

from chem_spectra.lib.converter.fid.base import FidBaseConverter
from chem_spectra.lib.converter.share import parse_params, parse_solvent

def search_brucker_processed(td):
    try:
        pdata_dir = find_dir(td, 'pdata')
        target_dir = get_processed_data(pdata_dir)
        return target_dir
    except:     # noqa: E722
        return False

def find_dir(path, name):
    for root, dirs, _ in os.walk(path):
        if name in dirs:
            return os.path.join(root, name)
    return False

def get_processed_data(path):
    processed_dirs = False
    for root, dirs, files in os.walk(path):
        if (len(dirs) > 0):
            processed_dirs = []
            for dir_name in dirs:
                dir_path = os.path.join(root, dir_name)
                processed_dirs.append(dir_path)
    return processed_dirs

class FidHasBruckerProcessed:
    def __init__(self, target_dir, params=False, fname=''):
        self.params = parse_params(params)
        self.data = self.__read(target_dir, fname)

    def __read(self, target_dir, fname):
        processed_dirs = search_brucker_processed(target_dir)
        data =[]

        unprocessed_dic, unprocessed_data = ng.bruker.read(target_dir)
        unprocessed_dic = self.__process_dic(unprocessed_dic, unprocessed_data, fname)
        unprocessed_data = self.__process_raw_data(unprocessed_dic, unprocessed_data)

        unprocessed_fid_conv = FidBaseConverter(dic=unprocessed_dic, data=unprocessed_data, params=self.params, fname=fname)
        data.append(unprocessed_fid_conv)

        for dir in processed_dirs:
            try:
                processed_dic, processed_data = ng.bruker.read_pdata(dir)
            except (OSError, FileNotFoundError):
                # skip silently if no binaries like 1r/1i
                continue

            processed_dic = self.__process_dic(processed_dic, processed_data, fname)
            
            try:
              processed_data = ng.process.proc_bl.baseline_corrector(processed_data, wd=20)  # baseline correction     # noqa: E501
            except:
              pass
            processed_data = ng.proc_base.di(processed_data)                # discard the imaginaries

            processed_fid_conv = FidBaseConverter(dic=processed_dic, data=processed_data, params=self.params, fname=fname)
            data.append(processed_fid_conv)
        return data

    def __process_dic(self, dic, data, fname):
        processed_dic = dic
        udic = ng.bruker.guess_udic(processed_dic, data).get(0) or {}
        # process dic
        processed_dic['.OBSERVENUCLEUS'] = '^{}'.format(udic.get('label'))
        processed_dic['.OBSERVEFREQUENCY'] = [udic.get('obs')]
        processed_dic['.SOLVENTNAME'] = processed_dic.get('acqus', {}).get('SOLVENT')
        processed_dic['.SHIFTREFERENCE'] = processed_dic.get('acqus', {}).get('SOLVENT')
        processed_dic['.PULSESEQUENCE'] = processed_dic.get('acqus', {}).get('PULPROG')
        offset = (float(processed_dic['acqus']['SW']) / 2) - (float(processed_dic['acqus']['O1']) / float(processed_dic['acqus']['BF1']))     # noqa: E501
        pt_head = float(processed_dic['acqus']['SW']) - offset
        pt_tail = -offset
        processed_dic['$OFFSET'] = [offset]
        processed_dic['FIRSTX'] = [pt_head]
        processed_dic['LASTX'] = [pt_tail]
        processed_dic['XUNITS'] = ['PPM']
        processed_dic['YUNITS'] = ['ARBITRARY']
        processed_dic['TITLE'] = ['FID {}'.format('.'.join(fname.split('.')[:-1]))]
        processed_dic['$CSSOLVENTX']     = [f'{offset:.6f}']
        processed_dic['$CSSOLVENTVALUE'] = ['0.000000']
        processed_dic['$CSSOLVENTNAME']  = ['AUTO-OFFSET']
        return processed_dic

    def __process_raw_data(self, dic, data):
        processed_data = data
        num_pts = processed_data.shape[-1]
        # process data (i.e. ys)
        processed_data = ng.bruker.remove_digital_filter(dic, processed_data)  # remove the digital filter   # noqa: E501
        processed_data = ng.proc_base.zf_size(processed_data, num_pts)    # zero fill to 32768 points   # noqa: E501
        processed_data = ng.proc_base.fft(processed_data)               # Fourier transform
        # p0, p1 = ng.process.proc_autophase.manual_ps(data)
        # data = ng.proc_base.ps(data, p0=-60, p1=200)
        if 'dept' not in dic['.PULSESEQUENCE']:
            data_am = ng.process.proc_autophase.autops(processed_data, 'acme')  # phase correction    # noqa: E501
            data_pm = ng.process.proc_autophase.autops(processed_data, 'peak_minima')  # phase correction     # noqa: E501
            processed_data = data_am if (data_am.min() > data_pm.min()) else data_pm
        else:
            processed_data = ng.proc_base.ps(processed_data, p0=0, p1=210)
            data_pm = ng.process.proc_autophase.autops(processed_data, 'peak_minima')  # phase correction     # noqa: E501
            processed_data = data_pm
        processed_data = ng.process.proc_bl.baseline_corrector(processed_data, wd=20)  # baseline correction     # noqa: E501
        processed_data = ng.proc_base.di(processed_data)                # discard the imaginaries
        processed_data = ng.proc_base.rev(processed_data)               # reverse the data
        return processed_data
 
