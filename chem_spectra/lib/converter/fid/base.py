import nmrglue as ng
import numpy as np

from chem_spectra.lib.converter.share import parse_params


class FidBaseConverter:
    def __init__(self, target_dir, params=False, fname=''):
        self.params = parse_params(params)
        self.dic, self.data = self.__read(target_dir, fname)
        self.datatypes = ['NMR SPECTRUM']
        self.datatype = 'NMR SPECTRUM'
        self.title = self.dic.get('TITLE', [''])[0]
        self.typ = 'NMR'
        self.fname = '.'.join(params.get('fname').split('.')[:-1])
        self.is_em_wave = self.__is_em_wave()
        self.is_ir = self.__is_ir()
        self.ncl = self.__ncl()
        self.simu_peaks = self.__read_simu_peaks()
        self.solv_peaks = []
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
        return self.typ in ['INFRARED', 'RAMAN']

    def __is_ir(self):
        return self.typ in ['INFRARED']

    def __ncl(self):
        try:
            ncls = self.dic['NUC1']
            if '^1H' in ncls:
                return '1H'
            elif '13C' in ncls:
                return '13C'
            elif '19F' in ncls:
                return '19F'
        except: # noqa
            pass
        return '1H'

    def __read_simu_peaks(self):
        target = self.dic.get('$CSSIMULATIONPEAKS', [])
        if target:
            target = [float(t) for t in target[0].split('\n')]
            return sorted(target)
        return []

    def __read_solvent(self):
        if self.ncl == '13C':
            ref_name = (
                self.params['ref_name'] or
                self.dic.get('$CSSOLVENTNAME', [''])[0]
            )
            # if ref_name and ref_name != '- - -':
            if ref_name: # skip when the solvent is exist.
                return
            orig_solv = (
                self.dic.get('.SOLVENTNAME', [''])[0] + \
                self.dic.get('.SHIFTREFERENCE', [''])[0]
            ).lower()

            if 'acetone' in orig_solv:
                self.dic['$CSSOLVENTNAME'] = ['Acetone-d6 (sep)']
                self.dic['$CSSOLVENTVALUE'] = ['29.920']
                self.dic['$CSSOLVENTX'] = ['0']
                self.solv_peaks = [(27.0, 33.0), (203.7, 209.7)]
