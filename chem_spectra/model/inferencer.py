import io
import requests
import numpy as np
import json
from flask import current_app

from chem_spectra.model.molecule import MoleculeModel
from chem_spectra.lib.data_pipeline.infrared import InfraredModel

hdr_nsdb = {
    'Content-Type': 'application/json'
}


class InferencerModel:
    def __init__(
        self,
        molfile=False, layout=False, peaks=False, shift=False, spectrum=False
    ):
        self.molfile = molfile.core
        self.layout = layout
        self.peaks = peaks
        self.shift = shift
        self.spectrum = spectrum


    @classmethod
    def predict_nmr(
        cls,
        molfile=False, layout=False, peaks=False, shift=False
    ):
        instance = cls(
            molfile=molfile,
            layout=layout,
            peaks=peaks,
            shift=shift
        )
        return instance.__predict_nmr()


    def __predict_nmr(self):
        peak_xs = self.__extract_x()
        solvent = self.shift.get('ref', {}) .get('nsdb')

        if self.layout == '1H':
            typ = 'nmr;1H;1d'
            data = self.__build_data(typ, peak_xs, solvent)
            rsp = requests.post(
                current_app.config['URL_NSHIFTDB'],
                headers=hdr_nsdb,
                json=data,
            )
            return rsp
        elif self.layout == '13C':
            typ = 'nmr;13C;1d'
            data = self.__build_data(typ, peak_xs, solvent)
            rsp = requests.post(
                current_app.config['URL_NSHIFTDB'],
                headers=hdr_nsdb,
                json=data,
            )
            return rsp

        return False


    def __extract_x(self):
        total = []
        for p in self.peaks:
            total.append(str(p['x']))

        return ';'.join(total)


    def __build_data(self, typ, peak_xs, solvent):
        return {
            'inputs':[
                {
                    'id': 1,
                    'type': typ,
                    'shifts': peak_xs,
                    'solvent': solvent,
                },
            ],
            'moltxt': self.molfile
        }


    @classmethod
    def predict_ir(cls, molfile=False, spectrum=False):
        instance = cls(
            molfile=molfile,
            spectrum=spectrum
        )
        return instance.__predict_ir()


    def __predict_ir(self):
        mm = MoleculeModel(self.molfile)
        fgs = { 'fgs': json.dumps(mm.fgs()) }

        im = InfraredModel(self.spectrum)
        xs, ys = im.standarize()

        buf = io.BytesIO()
        np.savez(buf, ys=ys)
        file = buf.getvalue()
        files = { 'file': (file) }

        rsp = requests.post(
            current_app.config['URL_DEEPIR'],
            files=files,
            data=fgs,
        )
        return rsp
