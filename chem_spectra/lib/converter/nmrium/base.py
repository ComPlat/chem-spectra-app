import nmrglue as ng
import json
import numpy as np
from chem_spectra.lib.converter.datatable import DatatableModel
from chem_spectra.lib.shared.calc import (  # noqa: E402
    calc_mpy_center
)

def coupling_string(js):
    if len(js) == 0:
        return ''
    return ', ' + ' '.join([str(j) for j in js])

class NMRiumDataConverter:
    def __init__(self, file=None):
        self.params = {'integration':{}, 'multiplicity':{}, 'ref_name':'', 'ref_value':'', 'select_x':0}
        self.file = file
        self.fname = ''
        self.is_em_wave = False
        self.non_nmr = False
        self.data = None

        spectrum_data = self.__read_file()
        self.edit_peaks = []
        self.auto_peaks = []

        if spectrum_data is not None:
            self.data = self.__parsing_xy_values(spectrum_data)
            self.__read_info(spectrum_data)

            self.boundary = self.__find_boundary()
            self.datatable = self.__set_datatable()

            self.mpy_itg_table, self.mpy_pks_table = self.__read_multiplicity(spectrum_data)
            self.edit_peaks = self.__read_peaks(spectrum_data)
        
        self.itg_table = []
        self.simu_peaks = []
        self.is_cyclic_volta = False
        self.typ = ''
        self.threshold = 1.0
        
        

    def __read_file(self):
        if (self.file is None):
            return None
        
        try:
            rawData = json.loads(self.file.core)
        except ValueError as e:
            return None
        
        parsedData = self.__parsing_spectra(rawData)

        return parsedData


    def __parsing_spectra(self, jsonData=None):
        if jsonData is None:
            return None
        
        spectra = jsonData['spectra']
        displaying_spectra = self.__find_displaying_spectra(spectra)
        numberOfSpectrum = len(displaying_spectra)

        if numberOfSpectrum <= 0:
            return None

        correlations = jsonData['correlations']
        displaying_spectrum_id = self.__find_displaying_spectrum_id(correlations)

        displayingSpectrum = None

        if (displaying_spectrum_id == '' and len(displaying_spectra)):
            displayingSpectrum = displaying_spectra[0]
        else:
            for spectrum in displaying_spectra:
                if spectrum['id'] == displaying_spectrum_id:
                    displayingSpectrum = spectrum
                    break

        return displayingSpectrum

    def __find_displaying_spectra(self, spectra=None):
        if spectra is None:
            return None

        filtered_spectra = filter(self.__check_displaying_spectrum, spectra)
        filtered_spectra = list(filtered_spectra)
        
        return filtered_spectra

    def __find_displaying_spectrum_id(self, correlations=None):
        if correlations is None:
            return ''
        
        try:
            values = correlations['values']
            if len(values) > 0:
                first_value = values[0]
                link = first_value['link']
                if len(link) > 0:
                    first_link = link[0]
                    experimentID = first_link['experimentID']
                    return experimentID
            return ''
        except Exception as e:
            return ''

    def __check_displaying_spectrum(self, spectrum):
        try:
            return spectrum['info']['isFid'] == False
        except Exception as e:
            return False

    def __parsing_xy_values(self, spectrumData):
        if spectrumData is None:
            return None
        x_values, y_values = self.__read_xy_values(spectrumData)

        self.xs = np.array(x_values)
        self.ys = np.array(y_values)

        self.first_x = x_values[0]
        self.last_x = x_values[len(x_values)-1]
        self.data_format = '(XY..XY)'

        return {"x":x_values, "y":y_values}


    def __read_info(self, spectrumData):
        if spectrumData is None:
            return None
        self.fname = self.file.name
        
        self.label = {"x": "ppm", "y": "intensity"}
        self.factor = {"x": 1.0, "y": 1.0}

        dic_info = spectrumData['info']
        self.datatype = dic_info['type']

        dic_meta = spectrumData['meta']
        self.dic = dic_meta

    def __read_xy_values(self, spectrum):
        if spectrum is None:
            return None, None

        data = spectrum['data']
        x_values = data['x']
        y_values = data['re']

        return x_values, y_values

    def __find_boundary(self):
        return {
            'x': {
                'max': self.xs.max(),
                'min': self.xs.min(),
            },
            'y': {
                'max': self.ys.max(),
                'min': self.ys.min(),
            },
        }

    def __set_datatable(self):
        y_factor = self.factor and self.factor['y']
        y_factor = y_factor or 1.0
        if (self.data_format and self.data_format == '(XY..XY)'):
            return DatatableModel().encode(
                self.ys,
                y_factor,
                self.xs,
                True
            )
        return DatatableModel().encode(
            self.ys,
            y_factor
        )

    def __read_multiplicity(self, spectrumData):
        ascii_uppercase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        dic_ranges = spectrumData['ranges']
        arr_values = dic_ranges['values']
        mpy_itg_table = []
        arr_multiplicites = []

        mpy_pks_table = []
        for rangeIdx, rangeVal in enumerate(arr_values):
            signalData = self.__read_signal(rangeVal['signals'])
            rangeData = self.__read_range_data(rangeVal)
            if ((signalData is not None) and (rangeData is not None)):
                for signal in signalData:
                    mpy_item = {
                        'mpyType': signal['mpyType'], 
                        'xExtent': rangeData['xExtent'],
                        'yExtent': rangeData['yExtent'],
                        'peaks': signal['peaks'],
                        'area': signal['area'],
                        'js': ''}
                    arr_multiplicites.append(mpy_item)

        for idx, mpy in enumerate(arr_multiplicites):
            mpy_itg_table.extend([
                '({}, {}, {}, {}, {}, {}, {}, {}{})\n'.format(
                    idx + 1,
                    mpy['xExtent']['xL'],
                    mpy['xExtent']['xU'] ,
                    calc_mpy_center(mpy['peaks'], 0, mpy['mpyType']),   # noqa: E501
                    float(mpy['area']) * 1.0,
                    idx + 1,
                    mpy['mpyType'],
                    ascii_uppercase[idx],
                    coupling_string(mpy['js']),
                ),
            ])

            for p in mpy['peaks']:
                mpy_pks_table.extend([
                    '({}, {}, {})\n'.format(
                                idx+1,
                                p['x'],
                                p['y'],
                            ),
                ])
        return mpy_itg_table, mpy_pks_table

    def __read_signal(self, signals):
        if len(signals) <= 0:
            return None

        parsedData = []
        for signalItem in signals:
            peaks = signalItem['peaks']
            mpyType = signalItem['multiplicity']
            area = signalItem['integration']
            parsedData.append({"mpyType": mpyType, "peaks": peaks, "area": area})

        return parsedData

    def __read_range_data(self, rangeValue):
        if rangeValue is None:
            return None

        xExtent = {"xL": rangeValue['from'], "xU": rangeValue['to']}
        yExtent = {"yL": 0.0, "yU": 0.0}
        return {"xExtent": xExtent, "yExtent": yExtent}
    
    def __read_peaks(self, spectrumData):
        dic_ranges = spectrumData['peaks']
        arr_values = dic_ranges['values']
        peaks_x = []
        peaks_y = []
        for peak_value in arr_values:
            peaks_x.append(peak_value['x'])
            peaks_y.append(peak_value['y'])
        return { 'x': peaks_x, 'y': peaks_y }
        
        

