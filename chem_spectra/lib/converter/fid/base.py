import nmrglue as ng
import numpy as np

from chem_spectra.lib.converter.share import parse_params, parse_solvent


class FidBaseConverter:
    def __init__(self, target_dir, params=False, fname=''):
        self.params = parse_params(params)
        self.dic, self.data = self.__read(target_dir, fname)
        self.datatypes = ['NMR SPECTRUM']
        self.datatype = 'NMR SPECTRUM'
        self.dataclass = None
        self.data_format = None
        self.title = self.dic.get('TITLE', [''])[0]
        self.typ = 'NMR'
        self.fname = '.'.join(params.get('fname').split('.')[:-1])
        self.is_em_wave = self.__is_em_wave()
        self.is_ir = self.__is_ir()
        self.is_tga = self.__is_tga()
        self.is_uv_vis = self.__is_uv_vis()
        self.is_hplc_uv_vis = self.__is_hplc_uv_vis()
        self.is_xrd = self.__is_xrd()
        self.is_cyclic_volta = self.__is_cyclic_volta()
        self.non_nmr = self.__non_nmr()
        self.ncl = self.__ncl()
        self.simu_peaks = self.__read_simu_peaks()
        self.solv_peaks = []
        self.is_dept = self.__is_dept()
        self.__read_solvent()

    def __read(self, target_dir, fname):
        dic, data = ng.bruker.read(target_dir)
        udic = ng.bruker.guess_udic(dic, data).get(0) or {}
        # process dic
        dic['.OBSERVENUCLEUS'] = '^{}'.format(udic.get('label'))
        dic['.OBSERVEFREQUENCY'] = [udic.get('obs')]
        dic['.SOLVENTNAME'] = dic.get('acqus', {}).get('SOLVENT')
        dic['.SHIFTREFERENCE'] = dic.get('acqus', {}).get('SOLVENT')
        dic['.PULSESEQUENCE'] = dic.get('acqus', {}).get('PULPROG')
        num_pts = data.shape[-1]
        offset = (float(dic['acqus']['SW']) / 2) - (float(dic['acqus']['O1']) / float(dic['acqus']['BF1']))
        pt_head = float(dic['acqus']['SW']) - offset
        pt_tail = -offset
        dic['$OFFSET'] = [offset]
        dic['FIRSTX'] = [pt_head]
        dic['LASTX'] = [pt_tail]
        dic['XUNITS'] = ['PPM']
        dic['YUNITS'] = ['ARBITRARY']
        dic['TITLE'] = ['FID {}'.format('.'.join(fname.split('.')[:-1]))]

        # process data (i.e. ys)
        data = ng.bruker.remove_digital_filter(dic, data) # remove the digital filter
        data = ng.proc_base.zf_size(data, num_pts)    # zero fill to 32768 points
        data = ng.proc_base.fft(data)               # Fourier transform
        # p0, p1 = ng.process.proc_autophase.manual_ps(data)
        # data = ng.proc_base.ps(data, p0=-60, p1=200)
        if 'dept' not in dic['.PULSESEQUENCE']:
            data_am = ng.process.proc_autophase.autops(data, 'acme') # phase correction
            data_pm = ng.process.proc_autophase.autops(data, 'peak_minima') # phase correction
            data = data_am if (data_am.min() > data_pm.min()) else data_pm
        else:
            data = ng.proc_base.ps(data, p0=0, p1=210)
            data_pm = ng.process.proc_autophase.autops(data, 'peak_minima') # phase correction
            data = data_pm
        data = ng.process.proc_bl.baseline_corrector(data, wd=20) # baseline correction
        data = ng.proc_base.di(data)                # discard the imaginaries
        data = ng.proc_base.rev(data)               # reverse the data
        return dic, data

    def __is_em_wave(self):
        return self.typ in ['INFRARED', 'RAMAN', 'UVVIS']

    def __non_nmr(self):
        return self.typ in ['INFRARED', 'RAMAN', 'UVVIS', 'HPLC UVVIS', 'THERMOGRAVIMETRIC ANALYSIS', 'MS', 'X-RAY DIFFRACTION', 'CYCLIC VOLTAMMETRY']

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

    def __ncl(self):
        try:
            ncls = self.dic.get('NUC1') or self.dic.get('.OBSERVENUCLEUS')
            if '^1H' in ncls:
                return '1H'
            elif '13C' in ncls:
                return '13C'
            elif '19F' in ncls:
                return '19F'
            elif '31P' in ncls:
                return '31P'
            elif '15N' in ncls:
                return '15N'
            elif '29Si' in ncls:
                return '29Si'
        except: # noqa
            pass
        return '13C'

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

        try: # TBD
            for p in (self.dic.get('.PULSESEQUENCE', []) + self.dic.get('.PULSE SEQUENCE', [])):
                if 'dept' in p:
                    return True
        except:
            pass

        return False
