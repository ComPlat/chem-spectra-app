import netCDF4

from chem_spectra.lib.converter.share import parse_params


class CdfBaseConverter:
    def __init__(self, path, params=False):
        self.params = parse_params(params)
        self.dic, self.data = self.__read(path)
        self.datatypes = ['MASS SPECTRUM']
        self.datatype = 'MASS SPECTRUM'
        self.title = self.dic.get('TITLE', [''])[0]
        self.typ = 'MS'

    def __read(self, path):
        cdf = netCDF4.Dataset(path)
        data = {
            'XS': cdf.variables['mass_values'],
            'YS': cdf.variables['intensity_values'],
            'SCANINDEXS': cdf.variables['scan_index'],
        }
        dic = {
            'TITLE': cdf.source_file_reference,
            'SPECTROMETER TYPE':cdf.experiment_type,
            'INLET': cdf.test_ms_inlet,
            'IONIZATION MODE': cdf.test_ionization_mode,
        }
        return dic, data
