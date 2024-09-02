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
        tic_positive_file_path = os.path.join(target_dir, 'TIC_PLUS.csv')
        tic_postive_data = self.__read_tic(tic_positive_file_path)

        tic_negative_file_path = os.path.join(target_dir, 'TIC_MINUS.csv')
        tic_negative_data = self.__read_tic(tic_negative_file_path, True)

        spectra_file_path = os.path.join(target_dir, 'MZ_Spectra.csv')
        data_frame = pd.read_csv(spectra_file_path, index_col='time', header=0)
        grouped_df = data_frame.groupby('time').agg(list)
        data_dict = {time: {'mz': group['mz'], 'intensities': group['intensities']} for time, group in grouped_df.iterrows()}
        return [tic_postive_data, tic_negative_data, data_dict]

    def __read_tic(self, file_path, is_negative = False):
        data_frame = pd.read_csv(file_path, header=0)
        tic_postive_data = data_frame.to_dict(orient='list')
        return tic_postive_data
        

