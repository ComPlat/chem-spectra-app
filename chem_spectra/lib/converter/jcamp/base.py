import nmrglue as ng
import json

from chem_spectra.lib.converter.share import parse_params, parse_solvent
from chem_spectra.lib.converter.jcamp.data_parse import read_parsed_jdx_data
import os

data_type_json = os.path.join(os.path.dirname(__file__), 'data_type.json')


class JcampBaseConverter:
    def __init__(self, path, params=False):
        self.params = parse_params(params)
        self.dic, self.data = self.__read(path)
        self.datatypes = self.dic['DATATYPE']
        self.datatypes = [datatype.upper() for datatype in self.datatypes]
        self.datatype = self.__set_datatype()
        self.dataclasses = {}
        if 'DATACLASS' in self.dic:
            self.dataclasses = self.dic['DATACLASS']
        self.dataclass = self.__set_dataclass()
        self.data_format = self.__set_dataformat()
        self.title = self.dic.get('TITLE', [''])[0]
        self.typ = self.__typ()
        self.fname = self.params.get('fname')
        self.is_em_wave = self.__is_em_wave()
        self.is_ir = self.__is_ir()
        self.is_tga = self.__is_tga()
        self.is_uv_vis = self.__is_uv_vis()
        self.is_hplc_uv_vis = self.__is_hplc_uv_vis()
        self.is_xrd = self.__is_xrd()
        self.is_cyclic_volta = self.__is_cyclic_volta()
        self.is_sec = self.__is_sec()
        self.is_cds = self.__is_cds()
        self.is_aif = self.__is_aif()
        self.is_emissions = self.__is_emissions()
        self.is_dls_acf = self.__is_dls_acf()
        self.is_dls_intensity = self.__is_dls_intensity()
        self.non_nmr = self.__non_nmr()
        self.ncl = self.__ncl()
        self.simu_peaks = self.__read_simu_peaks()
        self.solv_peaks = []
        self.is_dept = self.__is_dept()
        self.__read_solvent()

    def __read(self, path):
        parsed_data = ng.jcampdx.read(path, show_all_data=True, read_err='ignore')
        return_dic, return_data = read_parsed_jdx_data(parsed_data)
        return return_dic, return_data

    def __set_datatype(self):
        dts = self.datatypes
        dt_dict = {
            'NMR': 'NMR SPECTRUM',
            'INFRARED': 'INFRARED SPECTRUM',
            'RAMAN': 'RAMAN SPECTRUM',
            'MS': 'MASS SPECTRUM',
            'HPLC UVVIS': 'HPLC UV/VIS SPECTRUM',
            'UVVIS': 'UV/VIS SPECTRUM',
        }

        with open(data_type_json, 'r') as mapping_file:
            data_type_mappings = json.load(mapping_file)["datatypes"]
        for key, values in data_type_mappings.items():
            values = [value.upper() for value in values]
            for dt in dts:
                if dt in values and key in dt_dict:
                    return dt_dict[key]
                elif dt in values and not key in dt_dict:
                    return key
        return ''

    def __typ(self):
        dt = self.datatype

        with open(data_type_json, 'r') as mapping_file:
            data_type_mappings = json.load(mapping_file)["datatypes"]

        for key, values in data_type_mappings.items():
            values = [value.upper() for value in values]
            if dt.upper() in values:
                return key
        return ''

    def __set_dataclass(self):
        data_class = self.dataclasses
        if 'XYPOINTS' in data_class:
            return 'XYPOINTS'
        elif 'XYDATA' in data_class:
            return 'XYDATA_OLD'
        return ''

    def __set_dataformat(self):
        try:
            return self.dic[self.dataclass][0].split('\n')[0]
        except: # noqa
            pass
        return '(X++(Y..Y))'

    def __is_em_wave(self):
        return self.typ in ['INFRARED', 'RAMAN', 'UVVIS']

    def __non_nmr(self):
        with open(data_type_json, 'r') as mapping_file:
            data_type_mappings = json.load(mapping_file).get("datatypes")
        dts = [dt for dt in data_type_mappings.keys() if dt != 'NMR']
        return self.typ in dts

    def __is_ir(self):
        return self.typ in ['INFRARED']

    def __is_tga(self):
        return self.typ in ['THERMOGRAVIMETRIC ANALYSIS']

    def __is_uv_vis(self):
        return self.typ in ['UVVIS']

    def __is_hplc_uv_vis(self):
        return self.typ in ['HPLC UVVIS']

    def __is_xrd(self):
        return self.typ in ['X-RAY DIFFRACTION']

    def __is_cyclic_volta(self):
        return self.typ in ['CYCLIC VOLTAMMETRY']

    def __is_sec(self):
        return self.typ in ['SIZE EXCLUSION CHROMATOGRAPHY']
    
    def __is_cds(self):
        return self.typ in ['CIRCULAR DICHROISM SPECTROSCOPY']

    def __is_aif(self):
        return self.typ in ['SORPTION-DESORPTION MEASUREMENT']
        
    def __is_emissions(self):
        return self.typ in ['Emissions']
    
    def __is_dls_acf(self):
        return self.typ in ['DLS ACF']

    def __is_dls_intensity(self):
        return self.typ in ['DLS intensity']

    def __ncl(self):
        try:
            ncls = self.dic['.OBSERVENUCLEUS']
            if '^1H' in ncls:
                return '1H'
            elif '^13C' in ncls:
                return '13C'
            elif '^19F' in ncls:
                return '19F'
            elif '31P' in ncls:
                return '31P'
            elif '15N' in ncls:
                return '15N'
            elif '29Si' in ncls:
                return '29Si'
        except: # noqa
            pass
        return ''

    def __read_simu_peaks(self):
        target = self.dic.get('$CSSIMULATIONPEAKS', [])
        if target:
            target = [float(t) for t in target[0].split('\n')]
            return sorted(target)
        return []

    def __read_solvent(self):
        parse_solvent(self)

    def __is_dept(self):
        if not self.ncl == '13C':
            return False

        for p in (self.dic.get('.PULSESEQUENCE', []) + self.dic.get('.PULSE SEQUENCE', [])):
            if 'dept' in p:
                return True

        return False
