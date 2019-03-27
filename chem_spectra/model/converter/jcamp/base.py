import nmrglue as ng
import numpy as np
from scipy import signal

from chem_spectra.model.converter.encoder import encode_datatable


class JcampBaseConverter():
    def __init__(self, path, params=False):
        self.params = self.__set_params(params)
        self.dic, self.data = self.__read(path)
        self.datatypes = self.dic['DATATYPE']
        self.datatype = self.__set_datatype()
        self.title = self.dic.get('TITLE', [''])[0]
        self.typ = self.__typ()


    def __set_params(self, params):
        if not params:
            return {
                'select_x': None,
                'ref_name': None,
                'ref_value': None,
                'peaks_str': None,
                'delta': 0.0,
                'mass': 0,
                'scan': None,
                'thres': None,
            }

        select_x = params.get('select_x', None)
        ref_name = params.get('ref_name', None)
        ref_value = params.get('ref_value', None)
        peaks_str = params.get('peaks_str', None)
        delta = 0.0
        mass = params.get('mass', 0)
        scan = params.get('scan', None)
        thres = params.get('thres', None)

        try:
            if select_x and float(select_x) != 0.0 and ref_name != '- - -' :
                delta = float(ref_value) - float(select_x)
        except:
            pass

        return {
            'select_x': select_x,
            'ref_name': ref_name,
            'ref_value': ref_value,
            'peaks_str': peaks_str,
            'delta': delta,
            'mass': mass,
            'scan': scan,
            'thres': thres,
        }


    def __read(self, path):
        return ng.jcampdx.read(path, show_all_data=True, read_err='ignore')


    def __set_datatype(self):
        dts = self.datatypes
        if 'NMR SPECTRUM' in dts:
            return 'NMR SPECTRUM'
        elif 'NMRSPECTRUM' in dts: # MNova
            return 'NMR SPECTRUM'
        elif 'INFRARED SPECTRUM' in dts:
            return 'INFRARED SPECTRUM'
        elif 'MASS SPECTRUM' in dts:
            return 'MASS SPECTRUM'
        return ''


    def __typ(self):
        dt = self.datatype
        if 'NMR SPECTRUM' == dt:
            return 'NMR'
        elif 'NMRSPECTRUM' == dt: # MNova
            return 'NMR'
        elif 'INFRARED SPECTRUM' == dt:
            return 'INFRARED' # TBD
        elif 'MASS SPECTRUM' == dt:
            return 'MS'
        return ''
