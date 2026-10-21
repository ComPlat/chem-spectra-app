import io
import requests
import numpy as np
import json
from flask import current_app
import logging

from chem_spectra.lib.data_pipeline.infrared import InfraredLib
from chem_spectra.lib.converter.jcamp.data_parse import UnparsableJcampData
from chem_spectra.lib.converter.share import UnconvertibleSpectrum
from chem_spectra.lib.chem.artist import ArtistLib

hdr_nsdb = {
    'Content-Type': 'application/json'
}

NMR_LAYOUTS = {'1H': 'nmr;1H;1d', '13C': 'nmr;13C;1d', '19F': 'nmr;19F;1d'}
DEFAULT_PREDICT_TIMEOUT = 120  # seconds

# `outline.code` values. The editor's forecast panel shows a fixed message per
# code and never shows `text`: 400 is "Peak & Element count mismatch", 500 and
# above "Service is not available", anything else "Unknown state". So 400 stays
# reserved for what it already meant, and `text` carries the reason for any
# client that shows it.
UNUSABLE_REQUEST = 422      # nothing a predictor could be asked about
NOT_AVAILABLE = 503         # not configured here, or not reachable
BAD_ANSWER = 502            # it answered, with an error or a shape we can't use
TIMED_OUT = 504


def outline(code, text):
    return {'outline': {'code': code, 'text': text}}


class PredictorNotConfigured(Exception):
    def __init__(self, what):
        super().__init__(
            '{} prediction is not configured on this server.'.format(what)
        )


def _predictor_url(key, what):
    url = current_app.config.get(key)
    if not url:
        raise PredictorNotConfigured(what)
    return url


def _timeout():
    return current_app.config.get('PREDICT_TIMEOUT', DEFAULT_PREDICT_TIMEOUT)


def _request_failed(err, service):
    """The outline for a request to `service` that did not get an answer."""
    if isinstance(err, requests.Timeout):
        return outline(TIMED_OUT, '{} did not answer within {} s.'.format(
            service, _timeout()))
    if isinstance(err, requests.ConnectionError):
        return outline(NOT_AVAILABLE, 'No Server available! Please try it later.')
    return outline(BAD_ANSWER, '{} could not be asked: {}'.format(service, err))


class InferencerModel:
    def __init__(
        self,
        mm=False, tm=False, layout=False, peaks=False, shift=False, spectrum=False      # noqa
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
            shift={'ref': {'label': False, 'name': '- - -', 'value': 0}, 'peak': False, 'enable': False}    # noqa: E501
        )
        try:
            rsp = instance.__predict_nmr(timeout=120)
            output = rsp.json()
            logger = logging.getLogger(__name__)
            logger.setLevel(logging.INFO)
            log_msg = 'smiles: {smi} - nmrshiftdb_result: {nmrshiftdb}'.format(smi=mm.smi, nmrshiftdb=output)   # noqa: E501
            logger.info(log_msg)
            simulations = sorted([shift['prediction'] for shift in output['result'][0]['shifts']])  # noqa: E501
            return simulations
        except Exception as error:
            logger = logging.getLogger(__name__)
            logger.setLevel(logging.ERROR)
            log_msg = 'smiles: {smi} - nmrshiftdb_result: {nmrshiftdb}'.format(smi=mm.smi, nmrshiftdb=error)    # noqa: E501
            logger.error(log_msg)
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

        if layout not in NMR_LAYOUTS:
            return outline(
                UNUSABLE_REQUEST,
                'NMR prediction covers {}, not {}.'.format(
                    ', '.join(NMR_LAYOUTS), layout),
            )

        instance = cls(
            mm=mm,
            layout=layout,
            peaks=peaks,
            shift=shift
        )
        try:
            rsp = instance.__predict_nmr()
        except PredictorNotConfigured as err:
            return outline(NOT_AVAILABLE, str(err))
        except requests.RequestException as err:
            logging.getLogger(__name__).error('NMRShiftDB: %s', err)
            return _request_failed(err, 'NMRShiftDB')
        if rsp.status_code >= 500:
            return outline(BAD_ANSWER, 'NMRShiftDB answered HTTP {}.'.format(
                rsp.status_code))
        try:
            output = rsp.json()
            shifts = output['result'][0]['shifts']
        except json.decoder.JSONDecodeError:
            # unchanged: an answer that is not JSON has always been reported
            # as the mismatch, which is what the editor shows for 400
            return outline(
                400, 'Peak & Element count mismatch! Please check peak-picking.')
        except (KeyError, IndexError, TypeError):
            return outline(
                BAD_ANSWER, 'NMRShiftDB answered in a shape this server does '
                'not understand.')
        output['result'][0]['svgs'] = ArtistLib.draw_nmr(
            mm=mm,
            layout=layout,
            predictions=shifts,
        )
        return {
            'outline': {
                'code': 200,
                'text': 'NMR prediction success.',  # noqa
            },
            'output': output,
        }

    def __predict_nmr(self, timeout=None):
        url = _predictor_url('URL_NSHIFTDB', 'NMR')
        peak_xs = self.__extract_x()
        solvent = self.shift.get('ref', {}) .get('nsdb')
        data = self.__build_data(NMR_LAYOUTS[self.layout], peak_xs, solvent)
        return requests.post(
            url,
            headers=hdr_nsdb,
            json=data,
            timeout=timeout or _timeout(),
            verify=False
        )

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
            "inputs": [
                {
                    "id": 1,
                    "type": typ,
                    "shifts": peak_xs,
                    "solvent": str(solvent),
                },
            ],
            "moltxt": self.mm.moltxt
        }

    @classmethod
    def predict_ir(cls, mm=False, spectrum=False):
        instance = cls(
            mm=mm,
            spectrum=spectrum
        )
        try:
            outcome = instance.__predict_ir()
        except PredictorNotConfigured as err:
            return outline(NOT_AVAILABLE, str(err))
        except UnconvertibleSpectrum as err:
            # A spectrum the converter refuses (no readable x range, ...).
            # This endpoint reports errors in the body, so keep that contract
            # and pass the reason on rather than answering 422.
            return outline(400, 'IR Spectrum error!\n{}'.format(err))
        except (TypeError, UnparsableJcampData):
            # UnparsableJcampData: the uploaded spectrum has no readable data
            # array, which is a bad upload rather than a server fault
            return outline(
                400, 'IR Spectrum error!\nPlease feedback to System Admins.')
        except requests.exceptions.JSONDecodeError:
            # before RequestException, which it also is
            return outline(
                BAD_ANSWER, 'The IR predictor answered with something that '
                'is not JSON.')
        except requests.RequestException as err:
            logging.getLogger(__name__).error('IR predictor: %s', err)
            return _request_failed(err, 'The IR predictor')
        try:
            fgs = outcome['output']['result'][0]['fgs']
        except (KeyError, IndexError, TypeError):
            return outline(
                BAD_ANSWER, 'The IR predictor answered in a shape this server '
                'does not understand.')
        outcome['output']['result'][0]['svgs'] = ArtistLib.draw_ir(
            mm=mm,
            layout='IR',
            predictions=fgs,
        )
        return outcome

    def __predict_ir(self):
        fgs = {'fgs': json.dumps(self.mm.fgs())}

        # the spectrum first: a problem with it is worth reporting even where
        # no predictor is configured
        im = InfraredLib(self.spectrum)
        xs, ys = im.standarize()
        url = _predictor_url('URL_DEEPIR', 'IR')

        buf = io.BytesIO()
        np.savez(buf, ys=ys)
        file = buf.getvalue()
        files = {'file': (file)}

        rsp = requests.post(
            url,
            files=files,
            data=fgs,
            timeout=_timeout(),
        )
        if rsp.status_code >= 500:
            raise requests.HTTPError(
                'it answered HTTP {}'.format(rsp.status_code))
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
