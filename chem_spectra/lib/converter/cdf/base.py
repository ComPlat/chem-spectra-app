import netCDF4

from chem_spectra.lib.converter.share import parse_params


class CdfBaseConverter:
    def __init__(self, path, params=False):
        self.params = parse_params(params)
        self.dic, self.data = self.__read(path)
        self.datatypes = ['MASS SPECTRUM']
        self.datatype = 'MASS SPECTRUM'
        self.dataclass = None
        self.data_format = None
        self.title = self.dic.get('TITLE', [''])[0]
        self.typ = 'MS'
        self.fname = '.'.join(params.get('fname').split('.')[:-1])

    def __read(self, path):
        cdf = netCDF4.Dataset(path)
        data = {
            'XS': cdf.variables['mass_values'],
            'YS': cdf.variables['intensity_values'],
            'SCANINDEXS': cdf.variables['scan_index'],
        }
        d = cdf.__dict__
        dic = {
            'TITLE': d.get('source_file_reference', ['']),
            'SPECTROMETER TYPE': d.get('experiment_type'),
            'INLET': d.get('test_ms_inlet'),
            'IONIZATION MODE': d.get('test_ionization_mode'),
        }
        return dic, data
