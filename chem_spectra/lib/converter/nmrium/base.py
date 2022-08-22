import nmrglue as ng
import json
import numpy as np
from chem_spectra.lib.converter.datatable import DatatableModel

class NMRiumDataConverter:
    def __init__(self, file):
        self.file = file
        self.data = self.__read_file()


        self.params = {'integration':{}, 'multiplicity':{}, 'ref_name':'', 'ref_value':'', 'select_x':0}
        self.itg_table = []
        self.mpy_itg_table = []
        self.mpy_pks_table = []
        self.simu_peaks = []
        self.is_cyclic_volta = False
        self.typ = ''
        self.threshold = 1.0
        self.edit_peaks = []
        self.auto_peaks = []
        

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
        numberOfSpectrum = len(spectra)

        if numberOfSpectrum <= 0:
            return None

        lastestSpectrum = spectra[numberOfSpectrum-1]

        self.__read_info(lastestSpectrum)

        x_values, y_values = self.__read_xy_values(lastestSpectrum)

        self.xs = np.array(x_values)
        self.ys = np.array(y_values)

        self.first_x = x_values[0]
        self.last_x = x_values[len(x_values)-1]
        self.data_format = '(XY..XY)'

        self.boundary = self.__find_boundary()
        self.datatable = self.__set_datatable()

        self.mpys = self.__read_multiplicity(lastestSpectrum)

        return {"x":x_values, "y":y_values}

    def __read_info(self, spectrumData):
        self.fname = self.file.name
        self.is_em_wave = False
        self.non_nmr = False
        self.label = {"x": "ppm", "y": "intensity"}
        self.factor = {"x": 1.0, "y": 1.0}

        dic_info = spectrumData['info']
        self.datatype = dic_info['type']

        dic_meta = spectrumData['meta']
        self.dic = dic_meta

    def __read_xy_values(self, spectrum):
        if spectrum is None:
            return None

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
        dic_ranges = spectrumData['ranges']
        arr_values = dic_ranges['values']
        mpy_itg_table = []
        for range_val in arr_values:
            signalData = self.__read_signal(range_val['signals'])
            rangeData = self.__read_range_data(range_val)
            if ((signalData is not None) and (rangeData is not None)):
                for signal in signalData:
                    mpy_item = {
                        'mpyType': signal['mpyType'], 
                        'xExtent': rangeData['xExtent'],
                        'yExtent': rangeData['yExtent'],
                        'peaks': signal['peaks'],
                        'area': signal['area']}
                    mpy_itg_table.append(mpy_item)

        return mpy_itg_table

    def __read_signal(self, signals):
        if len(signals) <= 0:
            return None

        parsedData = []
        for signalItem in signals:
            peaks = signalItem['peaks']
            mpyType = signalItem['multiplicity']
            area = signalItem['delta']
            parsedData.append({"mpyType": mpyType, "peaks": peaks, "area": area})

        return parsedData

    def __read_range_data(self, rangeValue):
        if rangeValue is None:
            return None

        xExtent = {"xL": rangeValue['from'], "xU": rangeValue['to']}
        yExtent = {"yL": 0.0, "yU": 0.0}
        return {"xExtent": xExtent, "yExtent": yExtent}
        

