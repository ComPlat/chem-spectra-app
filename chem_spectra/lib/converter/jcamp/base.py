import nmrglue as ng
import json
import logging
import re

from chem_spectra.lib.converter.share import parse_params, parse_solvent
from chem_spectra.lib.converter.jcamp.techniques import technique_for
from chem_spectra.lib.converter.jcamp.technique import UnconvertibleSpectrum
import os

data_type_json = os.path.join(os.path.dirname(__file__), 'data_type.json')

logger = logging.getLogger(__name__)

# Enough to reach the records; a JCAMP header is a few hundred bytes and the
# rest of the file is numbers. The same window chemotion_ELN reads.
HEADER_BYTES = 64 * 1024

# JCAMP-DX 4.24 (5.1): a label is compared with spaces, dashes, underscores
# and slashes removed, and without regard to case. Anywhere in the label, not
# only between its words -- so `##N-UM D/IM=` is `##NUMDIM=`. Normalising the
# whole label is both closer to the spec and simpler than matching each
# spelling.
LABEL_NOISE_RE = re.compile(r'[\s_/-]')

# `nD` is the JCAMP-DX 6 spelling; some vendors write the number instead.
ND_VALUE_RE = re.compile(r'^\s*([2-9]|n)\s*D\s+NMR', re.IGNORECASE)


def header_records(header_text):
    """Every `##LABEL= value` in the text, label normalised, in file order."""
    for line in header_text.split('\n'):
        if not line.startswith('##'):
            continue
        label, sep, value = line[2:].partition('=')
        if not sep:
            continue
        yield LABEL_NOISE_RE.sub('', label).upper(), value.strip()


def declared_dimensions(header_text):
    """How many dimensions the header declares, or None if it does not say.

    Read from the text rather than from nmrglue's result. Not because nmrglue
    drops the records -- it keeps them on a small file -- but because a real
    2D dataset is large and need not parse at all, and there is no reason to
    spend that read on a file that will be refused either way.

    An explicit `##NUM DIM=` decides. Failing that, the datatype is read: it
    says `nD` or `2D` without always saying which n, so two is the least it
    can mean.
    """
    from_datatype = None
    for label, value in header_records(header_text):
        if label == 'NUMDIM':
            try:
                return int(value.split()[0])
            except (ValueError, IndexError):
                continue
        if label == 'DATATYPE' and from_datatype is None:
            match = ND_VALUE_RE.match(value)
            if match:
                head = match.group(1)
                from_datatype = 2 if head.lower() == 'n' else int(head)
    return from_datatype


def declares_nmr(header_text):
    """Whether the header says this is NMR at all.

    Only used to word the refusal -- an NMR file can be opened in NMRium,
    anything else cannot. The refusal itself does not depend on it: no
    technique here can hold a second dimension.
    """
    for label, value in header_records(header_text):
        if label == 'DATATYPE' and 'NMR' in value.upper():
            return True
        if label in ('.OBSERVENUCLEUS', 'OBSERVENUCLEUS'):
            return True
    return False


def read_header(path):
    """The first HEADER_BYTES of the file, as text.

    latin-1 because it maps every byte and so cannot raise: this runs before
    anything has established the file is even a JCAMP.
    """
    try:
        with open(path, 'rb') as handle:
            raw = handle.read(HEADER_BYTES).decode('latin-1')
    except (OSError, TypeError, ValueError):
        return ''
    # Line endings normalised so `^` finds a label whatever wrote the file.
    # A CR-only file -- classic Mac, and some instrument exports -- is one
    # long line to `re.MULTILINE`, so every record after the first was
    # invisible.
    return raw.replace('\r\n', '\n').replace('\r', '\n')

class JcampBaseConverter:
    def __init__(self, path, params=False):
        self.params = parse_params(params)
        self.__refuse_multi_dimensional(path)
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

    @staticmethod
    def __refuse_multi_dimensional(path):
        """A 2D dataset is not a curve, and this app has nowhere to put one.

        Read one row of it and you get a 1D spectrum with the acquisition
        time axis labelled as the spectrum's -- which is what happened, with
        a 200 and no indication that the answer meant nothing. Refusing is
        the whole fix: `xs`/`ys` cannot hold a matrix, so detecting and
        continuing would only move the failure.

        Checked before the file is parsed, in the one place every endpoint
        goes through, so /convert, /zip_jcamp, the save and refresh paths and
        each BagIt member are all covered by this single guard. A BagIt
        archive carrying one such member is refused whole, as an archive
        holding two bagits is -- processing part of an upload silently is the
        defect, not the remedy.
        """
        header = read_header(path)
        dimensions = declared_dimensions(header)
        if dimensions is None or dimensions <= 1:
            return
        if declares_nmr(header):
            raise UnconvertibleSpectrum(
                'this is a {}D NMR file. ChemSpectra reads one-dimensional '
                'spectra only; open it in NMRium instead'.format(dimensions)
            )
        raise UnconvertibleSpectrum(
            'this file declares {} dimensions. ChemSpectra reads '
            'one-dimensional spectra only'.format(dimensions)
        )

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
        # Nothing matched. Keep the file's own ##DATA TYPE= rather than
        # returning '', because the composer writes this value straight back
        # out (composer/technique.py) and 'DATATYPE' is suppressed from the
        # original-metadata dump (composer/base.py) -- so '' erased the only
        # record of what the file said it was. The spectrum still renders as a
        # generic curve either way; what is lost is the ability to reclassify
        # it later, which is exactly what happens when an under-specified
        # technique is added to data_type.json after the fact.
        return self.__unrecognised_datatype()

    # Blocks that JCAMP uses structurally, or that carry a derived table
    # rather than a measurement. None of them names the technique.
    #
    # Compared with spaces removed, because the same block is spelled both
    # ways in the wild: this app composes `NMRPEAKTABLE`, while
    # chemotion-converter-app emits `NMR PEAK TABLE`. The suffix rules cover
    # every per-technique variant of those two -- `INFRARED PEAK TABLE`,
    # `NMP PEAK ASSIGNMENTS` (its misspelling) and so on.
    #
    # `tests/lib/converter/jcamp/test_jcamp_datatype_classification.py
    # ::test_auxiliary_blocks_stay_unmapped` holds the authoritative list of
    # spellings, and pins that none of them is in data_type.json.
    AUXILIARY_DATATYPES = ('LINK', 'NMRFID', 'INFRAREDINTERFEROGRAM')
    AUXILIARY_SUFFIXES = ('PEAKTABLE', 'PEAKASSIGNMENTS')

    @classmethod
    def _is_auxiliary_datatype(cls, datatype):
        squashed = datatype.replace(' ', '')
        return (squashed in cls.AUXILIARY_DATATYPES
                or squashed.endswith(cls.AUXILIARY_SUFFIXES))

    def __unrecognised_datatype(self):
        for dt in self.datatypes:
            if self._is_auxiliary_datatype(dt):
                continue
            return dt
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

