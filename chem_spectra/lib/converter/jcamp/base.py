import nmrglue as ng
import json
import logging

from chem_spectra.lib.converter.share import parse_params, parse_solvent
from chem_spectra.lib.converter.jcamp.techniques import technique_for
import os

data_type_json = os.path.join(os.path.dirname(__file__), 'data_type.json')

logger = logging.getLogger(__name__)

class JcampBaseConverter:
    def __init__(self, path, params=False):
        self.params = parse_params(params)
        self.dic, self.data = self.__read(path)
        # A file with no ##DATA TYPE= at all raised KeyError straight out of
        # the request. An absent header is no more exceptional than an
        # unrecognised one, so it takes the same path.
        self.datatypes = self.dic.get('DATATYPE') or []
        self.datatypes = [datatype.upper() for datatype in self.datatypes]
        self.datatype = self.__set_datatype()
        self.dataclasses = {}
        if 'DATACLASS' in self.dic:
            self.dataclasses = self.dic['DATACLASS']
        self.dataclass = self.__set_dataclass()
        self.data_format = self.__set_dataformat()
        self.title = self.dic.get('TITLE', [''])[0]
        self.typ = self.__typ()
        self.fname = self.params.get('fname')
        if not self.typ:
            # a caller-supplied data_type_mapping REPLACES the built-in one,
            # so pointing at data_type.json would be useless advice there
            source = ('the data_type_mapping supplied with this request'
                      if self.params.get('user_data_type_mapping')
                      else 'data_type.json')
            logger.warning(
                'unrecognised ##DATA TYPE= %s in %r; processing it as a '
                'generic curve. Add it to %s if this app should handle it '
                'as a known technique.',
                self.datatypes, self.fname, source,
            )
        self.technique = technique_for(self.typ)
        self.ncl = self.__ncl()
        self.simu_peaks = self.__read_simu_peaks()
        self.solv_peaks = []
        self.__read_solvent()
        self.__read_user_data_type_mapping()

    def __read(self, path):
        return ng.jcampdx.read(path, show_all_data=True, read_err='ignore')
    
    def __read_user_data_type_mapping(self):
        user_dt_mapping = self.params.get('user_data_type_mapping')
        if user_dt_mapping == '' or user_dt_mapping is None:
            return ''
        else:
            return json.loads(user_dt_mapping)['datatypes']

    def __data_type_mappings(self):
        if self.params.get('user_data_type_mapping'):
            return self.__read_user_data_type_mapping()
        with open(data_type_json, 'r') as mapping_file:
            return json.load(mapping_file)['datatypes']

    def __set_datatype(self):
        dts = self.datatypes
        dt_dict = {
            'NMR': 'NMR SPECTRUM',
            'INFRARED': 'INFRARED SPECTRUM',
            'RAMAN': 'RAMAN SPECTRUM',
            'MS': 'MASS SPECTRUM',
            'HPLC UVVIS': 'HPLC UV/VIS SPECTRUM',
            'UVVIS': 'UV/VIS SPECTRUM',
        }

        data_type_mappings = self.__data_type_mappings()

        # The file's own block order decides, not the order data_type.json
        # happens to list its keys in. The first recognised ##DATA TYPE= is
        # the primary measurement; auxiliary blocks (NMR FID, peak tables)
        # are deliberately absent from the mapping so they are skipped here.
        for dt in dts:
            for key, values in data_type_mappings.items():
                if dt in [value.upper() for value in values]:
                    return dt_dict.get(key, key)
        return ''

    def __typ(self):
        dt = self.datatype

        data_type_mappings = self.__data_type_mappings()

        for key, values in data_type_mappings.items():
            values = [value.upper() for value in values]
            if dt.upper() in values:
                return key
        return ''

    # `non_nmr` is the one predicate that survives, and not as a shim.
    # JcampMSConverter copies it (converter/jcamp/ms.py) and carries no
    # descriptor, so BaseComposer._technique() reaches it through
    # `getattr(self.core, 'non_nmr', True)` -- that fallback is the MS path's
    # only route to a descriptor. It goes when MS is folded in as a
    # technique; see DEFERRED.md item 1. The other fifteen predicates are
    # gone: every decision they carried is a field on the descriptor.
    @property
    def non_nmr(self):
        return self.technique.key != 'NMR'

    def __set_dataclass(self):
        data_class = self.dataclasses
        if 'XYPOINTS' in data_class:
            return 'XYPOINTS'
        elif 'XYDATA' in data_class:
            return 'XYDATA_OLD'
        return ''

    def __set_dataformat(self):
        try:
            return self.dic[self.dataclass][0].split('\n')[0]
        except: # noqa
            pass
        return '(X++(Y..Y))'

    def __ncl(self):
        try:
            ncls = self.dic['.OBSERVENUCLEUS']
            if '^1H' in ncls:
                return '1H'
            elif '^13C' in ncls:
                return '13C'
            elif '^19F' in ncls:
                return '19F'
            elif '31P' in ncls:
                return '31P'
            elif '15N' in ncls:
                return '15N'
            elif '29Si' in ncls:
                return '29Si'
        except: # noqa
            pass
        return ''

    def __read_simu_peaks(self):
        target = self.dic.get('$CSSIMULATIONPEAKS', [])
        if target:
            target = [float(t) for t in target[0].split('\n')]
            return sorted(target)
        return []

    def __read_solvent(self):
        parse_solvent(self)

