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
        self.is_tga = self.__is_tga()
        self.is_uv_vis = self.__is_uv_vis()
        self.is_xrd = self.__is_xrd()
        self.non_nmr = self.__non_nmr()
        self.ncl = self.__ncl()
        self.simu_peaks = self.__read_simu_peaks()
        self.solv_peaks = []
        self.is_dept = self.__is_dept()
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
        elif 'UV/VIS SPECTRUM' in dts:
            return 'UV/VIS SPECTRUM'
        elif 'UV-VIS' in dts:
            return 'UV/VIS SPECTRUM'
        elif 'THERMOGRAVIMETRIC ANALYSIS' in dts:
            return 'THERMOGRAVIMETRIC ANALYSIS'
        elif 'X-RAY DIFFRACTION' in dts:
            return 'X-RAY DIFFRACTION'
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
        elif 'UV/VIS SPECTRUM' == dt:
            return 'UVVIS'
        elif 'UV-VIS' == dt:
            return 'UVVIS'
        elif 'THERMOGRAVIMETRIC ANALYSIS' == dt:
            return 'THERMOGRAVIMETRIC ANALYSIS'
        elif 'X-RAY DIFFRACTION' == dt:
            return 'X-RAY DIFFRACTION'
        return ''

    def __is_em_wave(self):
        return self.typ in ['INFRARED', 'RAMAN', 'UVVIS']

    def __non_nmr(self):
        return self.typ in ['INFRARED', 'RAMAN', 'UVVIS', 'THERMOGRAVIMETRIC ANALYSIS', 'MS', 'X-RAY DIFFRACTION']

    def __is_ir(self):
        return self.typ in ['INFRARED']

    def __is_tga(self):
        return self.typ in ['THERMOGRAVIMETRIC ANALYSIS']

    def __is_uv_vis(self):
        return self.typ in ['UVVIS']

    def __is_xrd(self):
        return self.typ in ['X-RAY DIFFRACTION']

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

