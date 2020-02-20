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

    def __read(self, target_dir, fname):
        dic, data = ng.bruker.read(target_dir)
        udic = ng.bruker.guess_udic(dic, data).get(0) or {}
        # process dic
        dic['.OBSERVENUCLEUS'] = '^{}'.format(udic.get('label'))
        dic['.OBSERVEFREQUENCY'] = [udic.get('obs')]
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
        data = ng.process.proc_autophase.autops(data, 'acme', p0=30, p1=-80) # phase correction
        data = ng.process.proc_bl.baseline_corrector(data, wd=20) # baseline correction
        data = ng.proc_base.di(data)                # discard the imaginaries
        data = ng.proc_base.rev(data)               # reverse the data
        return dic, data
