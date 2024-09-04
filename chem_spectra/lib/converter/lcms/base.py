import os
import pandas as pd
from chem_spectra.lib.converter.share import parse_params

class LCMSBaseConverter:
    def __init__(self, target_dir, params=False, fname=''):
        self.params = parse_params(params)
        self.typ = None
        self.fname = fname
        if target_dir is None:
            self.data = None
        else:
            self.data = self.__read(target_dir, fname)
  
    def __read(self, target_dir, fname):
        tic_postive_data = self.__read_tic(target_dir)
        tic_negative_data = self.__read_tic(target_dir, is_negative=True)

        uvvis_data = self.__read_uvvis(target_dir)

        spectra_postive_data = self.__read_spectra(target_dir)
        spectra_negative_data = self.__read_spectra(target_dir, is_negative=True)

        return [tic_postive_data, tic_negative_data, uvvis_data, spectra_postive_data, spectra_negative_data]

    def __read_tic(self, target_dir, is_negative = False):
        file_path = os.path.join(target_dir, 'TIC_PLUS.csv')
        if is_negative:
            file_path = os.path.join(target_dir, 'TIC_MINUS.csv')
        data_frame = pd.read_csv(file_path, header=0)
        tic_postive_data = data_frame.to_dict(orient='list')
        return tic_postive_data

    def __read_uvvis(self, target_dir):
        file_path = os.path.join(target_dir, 'LCMS.csv')
        data_frame = pd.read_csv(file_path, index_col='wavelength', header=0)
        grouped_df = data_frame.groupby('wavelength').agg(list)
        data_dict = {wavelength: {'RetentionTime': group['RetentionTime'], 'DetectorSignal': group['DetectorSignal']} for wavelength, group in grouped_df.iterrows()}
        return data_dict

    def __read_spectra(self, target_dir, is_negative = False):
        file_path = os.path.join(target_dir, 'MZ_PLUS_Spectra.csv')
        if is_negative:
            file_path = os.path.join(target_dir, 'MZ_MINUS_Spectra.csv')
        data_frame = pd.read_csv(file_path, index_col='time', header=0)
        grouped_df = data_frame.groupby('time').agg(list)
        data_dict = {time: {'mz': group['mz'], 'intensities': group['intensities']} for time, group in grouped_df.iterrows()}
        return data_dict
        

