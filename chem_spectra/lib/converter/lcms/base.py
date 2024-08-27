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
        spectra_file_path = os.path.join(target_dir, 'MZ_Spectra.csv')
        data_frame = pd.read_csv(spectra_file_path, index_col='time', header=0)
        grouped_df = data_frame.groupby('time').agg(list)
        grouped_dict = {time: {'mz': group['mz'], 'intensities': group['intensities']} for time, group in grouped_df.iterrows()}
        return grouped_dict
        

