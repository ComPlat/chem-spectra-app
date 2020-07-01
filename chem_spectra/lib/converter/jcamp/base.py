import nmrglue as ng

from chem_spectra.lib.converter.share import parse_params, parse_solvent


class JcampBaseConverter:
    def __init__(self, path, params=False):
        self.params = parse_params(params)
        self.dic, self.data = self.__read(path)
        self.datatypes = self.dic['DATATYPE']
        self.datatype = self.__set_datatype()
        self.title = self.dic.get('TITLE', [''])[0]
        self.typ = self.__typ()
        self.fname = self.params.get('fname')
        self.is_em_wave = self.__is_em_wave()
        self.is_ir = self.__is_ir()
        self.ncl = self.__ncl()
        self.simu_peaks = self.__read_simu_peaks()
        self.solv_peaks = []
        self.__read_solvent()

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
        elif 'RAMAN SPECTRUM' in dts:
            return 'RAMAN SPECTRUM'
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
        elif 'RAMAN SPECTRUM' == dt:
            return 'RAMAN'  # TBD
        elif 'MASS SPECTRUM' == dt:
            return 'MS'
        return ''

    def __is_em_wave(self):
        return self.typ in ['INFRARED', 'RAMAN']

    def __is_ir(self):
        return self.typ in ['INFRARED']

    def __ncl(self):
        try:
            ncls = self.dic['.OBSERVENUCLEUS']
            if '^1H' in ncls:
                return '1H'
            elif '^13C' in ncls:
                return '13C'
            elif '^19F' in ncls:
                return '19F'
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
