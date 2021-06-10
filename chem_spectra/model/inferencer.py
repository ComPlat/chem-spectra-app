import io
import requests
import numpy as np
import json
from flask import current_app
import logging

from chem_spectra.lib.data_pipeline.infrared import InfraredLib
from chem_spectra.lib.chem.artist import ArtistLib

hdr_nsdb = {
    'Content-Type': 'application/json'
}


class InferencerModel:
    def __init__(
        self,
        mm=False, tm=False, layout=False, peaks=False, shift=False, spectrum=False
    ):
        self.mm = mm
        self.tm = tm
        self.layout = layout
        self.peaks = peaks
        self.shift = shift
        self.spectrum = spectrum

    @classmethod
    def simulate_nmr(
        cls, mm=False, layout='13C'
    ):
        instance = cls(
            mm=mm,
            layout=layout,
            peaks=[{'y': 0, 'x': -9999}],
            shift={'ref': {'label': False, 'name': '- - -', 'value': 0}, 'peak': False, 'enable': False}
        )
        try:
            rsp = instance.__predict_nmr(timeout=20)
            output = rsp.json()
            logger = logging.getLogger(__name__)
            logger.setLevel(logging.INFO)
            log_msg = 'smiles: {smi} - nmrshiftdb_result: {nmrshiftdb}'.format(smi=mm.smi, nmrshiftdb=output)
            logger.info(log_msg)
            simulations = sorted([shift['prediction'] for shift in output['result'][0]['shifts']])
            return simulations
        except:
            return []

    @classmethod
    def predict_nmr(
        cls,
        mm=False, layout=False, peaks=False, shift=False
    ):
        if len(peaks) == 1 and peaks[0]['x'] == -1000:
            return {
                'outline': {
                    'code': 400,
                    'text': 'Peak & Element count mismatch! Please check peak-picking.',  # noqa
                }
            }

        instance = cls(
            mm=mm,
            layout=layout,
            peaks=peaks,
            shift=shift
        )
        try:
            rsp = instance.__predict_nmr()
            output = rsp.json()
            svgs = ArtistLib.draw_nmr(
                mm=mm,
                layout=layout,
                predictions=output['result'][0]['shifts'],
            )
            output['result'][0]['svgs'] = svgs
            return {
                'outline': {
                    'code': 200,
                    'text': 'NMR prediction success.',  # noqa
                },
                'output': output,
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

    def __predict_nmr(self, timeout=None):
        peak_xs = self.__extract_x()
        solvent = self.shift.get('ref', {}) .get('nsdb')

        if self.layout == '1H':
            typ = 'nmr;1H;1d'
            data = self.__build_data(typ, peak_xs, solvent)
            rsp = requests.post(
                current_app.config.get('URL_NSHIFTDB'),
                headers=hdr_nsdb,
                json=data,
                timeout=timeout,
            )
            return rsp
        elif self.layout == '13C':
            typ = 'nmr;13C;1d'
            data = self.__build_data(typ, peak_xs, solvent)
            rsp = requests.post(
                current_app.config.get('URL_NSHIFTDB'),
                headers=hdr_nsdb,
                json=data,
                timeout=timeout,
            )
            return rsp
        elif self.layout == '19F':
            typ = 'nmr;19F;1d'
            data = self.__build_data(typ, peak_xs, solvent)
            rsp = requests.post(
                current_app.config.get('URL_NSHIFTDB'),
                headers=hdr_nsdb,
                json=data,
                timeout=timeout,
            )
            return rsp
        return False

    def __extract_x(self):
        total = []
        deviation = 0.0001
        for p in self.peaks:
            target = str(p['x'])
            if target in total:
                target = str(p['x'] + deviation)
                deviation += 0.0001
            total.append(target)

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
            'moltxt': self.mm.moltxt
        }

    @classmethod
    def predict_ir(cls, mm=False, spectrum=False):
        instance = cls(
            mm=mm,
            spectrum=spectrum
        )
        try:
            outcome = instance.__predict_ir()
            svgs = ArtistLib.draw_ir(
                mm=mm,
                layout='IR',
                predictions=outcome['output']['result'][0]['fgs'],
            )
            outcome['output']['result'][0]['svgs'] = svgs
            return outcome
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
        fgs = {'fgs': json.dumps(self.mm.fgs())}

        im = InfraredLib(self.spectrum)
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
    def predict_ms(cls, mm=False, tm=False):
        instance = cls(
            mm=mm,
            tm=tm
        )
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
        bx, by, _, _, scan = self.tm.prism_peaks()
        return {
            'outline': {
                'code': 200,
            },
            'output': {
                'result': [
                    {
                        'type': 'ms',
                        'xs': bx,
                        'ys': by,
                        'thres': self.tm.core.thres,
                        'scan': scan,
                    }
                ],
            },
        }
