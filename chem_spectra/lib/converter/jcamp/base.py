import nmrglue as ng

from chem_spectra.lib.converter.share import parse_params


class JcampBaseConverter:
    def __init__(self, path, params=False):
        self.params = parse_params(params)
        self.dic, self.data = self.__read(path)
        self.datatypes = self.dic['DATATYPE']
        self.datatype = self.__set_datatype()
        self.title = self.dic.get('TITLE', [''])[0]
        self.typ = self.__typ()

    def __read(self, path):
        return ng.jcampdx.read(path, show_all_data=True, read_err='ignore')

    def __set_datatype(self):
        dts = self.datatypes
        if 'NMR SPECTRUM' in dts:
            return 'NMR SPECTRUM'
        elif 'NMRSPECTRUM' in dts:  # MNova
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
        elif 'NMRSPECTRUM' == dt:  # MNova
            return 'NMR'
        elif 'INFRARED SPECTRUM' == dt:
            return 'INFRARED'  # TBD
        elif 'MASS SPECTRUM' == dt:
            return 'MS'
        return ''
