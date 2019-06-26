import io
import requests
import numpy as np
import json
from flask import current_app

from chem_spectra.model.molecule import MoleculeModel
from chem_spectra.lib.data_pipeline.infrared import InfraredModel
from chem_spectra.model.transformer import TransformerModel as TraModel

hdr_nsdb = {
    'Content-Type': 'application/json'
}


class InferencerModel:
    def __init__(
        self,
        mm=False, layout=False, peaks=False, shift=False, spectrum=False
    ):
        self.moltxt = mm.moltxt
        self.layout = layout
        self.peaks = peaks
        self.shift = shift
        self.spectrum = spectrum

    @classmethod
    def predict_nmr(
        cls,
        mm=False, layout=False, peaks=False, shift=False
    ):
        instance = cls(
            mm=mm,
            layout=layout,
            peaks=peaks,
            shift=shift
        )
        try:
            rsp = instance.__predict_nmr()
            return {
                'outline': {
                    'code': 200,
                    'text': 'NMR prediction success.',  # noqa
                },
                'output': rsp.json(),
            }
        except json.decoder.JSONDecodeError:
            return {
                'outline': {
                    'code': 400,
                    'text': 'Peak & Element count mismatch! Please check peak-picking.',  # noqa
                }
            }
        except requests.ConnectionError:
            return {
                'outline': {
                    'code': 503,
                    'text': 'No Server available! Please try it later.',
                }
            }

    def __predict_nmr(self):
        peak_xs = self.__extract_x()
        solvent = self.shift.get('ref', {}) .get('nsdb')

        if self.layout == '1H':
            typ = 'nmr;1H;1d'
            data = self.__build_data(typ, peak_xs, solvent)
            rsp = requests.post(
                current_app.config.get('URL_NSHIFTDB'),
                headers=hdr_nsdb,
                json=data,
            )
            return rsp
        elif self.layout == '13C':
            typ = 'nmr;13C;1d'
            data = self.__build_data(typ, peak_xs, solvent)
            rsp = requests.post(
                current_app.config.get('URL_NSHIFTDB'),
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
            'inputs': [
                {
                    'id': 1,
                    'type': typ,
                    'shifts': peak_xs,
                    'solvent': solvent,
                },
            ],
            'moltxt': self.moltxt
        }

    @classmethod
    def predict_ir(cls, mm=False, spectrum=False):
        instance = cls(
            mm=mm,
            spectrum=spectrum
        )
        try:
            return instance.__predict_ir()
        except TypeError:
            return {
                'outline': {
                    'code': 400,
                    'text': 'IR Spectrum error!\nPlease feedback to System Admins.',  # noqa
                }
            }
        except requests.ConnectionError:
            return {
                'outline': {
                    'code': 503,
                    'text': 'No Server available! Please try it later.',
                }
            }

    def __predict_ir(self):
        mm = MoleculeModel(self.moltxt)
        fgs = {'fgs': json.dumps(mm.fgs())}

        im = InfraredModel(self.spectrum)
        xs, ys = im.standarize()

        buf = io.BytesIO()
        np.savez(buf, ys=ys)
        file = buf.getvalue()
        files = {'file': (file)}

        rsp = requests.post(
            current_app.config.get('URL_DEEPIR'),
            files=files,
            data=fgs,
        )
        return rsp.json()

    @classmethod
    def predict_ms(cls, mm=False, spectrum=False):
        instance = cls(
            mm=mm,
            spectrum=spectrum
        )
        return instance.__predict_ms()
        try:
            return instance.__predict_ms()
        except TypeError:
            return {
                'outline': {
                    'code': 400,
                    'text': 'MS Spectrum error!\nPlease feedback to System Admins.',  # noqa
                }
            }
        except requests.ConnectionError:
            return {
                'outline': {
                    'code': 503,
                    'text': 'No Server available! Please try it later.',
                }
            }

    def __predict_ms(self):
        cmpsr = TraModel(self.spectrum, {'ext': 'jdx'}).to_composer()
        bx, by, gx, gy, scan = cmpsr.prism_peaks()
        return json.dumps({
            'outline': {
                'code': 200,
            },
            'output': {
                'result': [
                    {
                        'type': 'ms',
                        'xs': bx,
                        'ys': by,
                        'thres': cmpsr.core.thres,
                        'scan': scan,
                    }
                ],
            },
        })
