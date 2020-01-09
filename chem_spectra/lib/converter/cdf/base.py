import netCDF4


class CdfBaseConverter:
    def __init__(self, path, params=False):
        self.params = self.__set_params(params)
        self.dic, self.data = self.__read(path)
        self.datatypes = ['MASS SPECTRUM']
        self.datatype = 'MASS SPECTRUM'
        self.title = self.dic.get('TITLE', [''])[0]
        self.typ = 'MS'

    def __set_params(self, params):
        if not params:
            return {
                'select_x': None,
                'ref_name': None,
                'ref_value': None,
                'peaks_str': None,
                'delta': 0.0,
                'mass': 0,
                'scan': None,
                'thres': None,
                'clear': False,
            }

        select_x = params.get('select_x', None)
        ref_name = params.get('ref_name', None)
        ref_value = params.get('ref_value', None)
        peaks_str = params.get('peaks_str', None)
        delta = 0.0
        mass = params.get('mass', 0)
        scan = params.get('scan', None)
        thres = params.get('thres', None)
        clear = params.get('clear', False)
        ext = params.get('ext', '')

        try:
            if select_x and float(select_x) != 0.0 and ref_name != '- - -':
                delta = float(ref_value) - float(select_x)
        except:  # noqa
            pass

        return {
            'select_x': select_x,
            'ref_name': ref_name,
            'ref_value': ref_value,
            'peaks_str': peaks_str,
            'delta': delta,
            'mass': mass,
            'scan': scan,
            'thres': thres,
            'clear': clear,
            'ext': ext,
        }

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
