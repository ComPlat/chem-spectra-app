import tempfile
from chem_spectra.lib.shared.misc import is_number
from chem_spectra.lib.shared.calc import (  # noqa: E402
    calc_mpy_center
)


def extrac_dic(core, key):
    query = core.dic.get(key, '')
    if type(query) is list:
        return query[0]
    return query


def calc_npoints(peaks):
    if peaks:
        return len(peaks['x'])
    return 0


def coupling_string(js):
    if len(js) == 0:
        return ''
    return ', ' + ' '.join([str(j) for j in js])


def is_metadata_to_be_ignored(keyword):
    return keyword in ['__comments', '_comments', 'FIRST', 'LAST', 'XYDATA_OLD', 'NTUPLES', 'PEAKASSIGNMENTS', 'XYDATA', '$CSSIMULATIONPEAKS', 'XFACTOR', 'YFACTOR', 'FIRSTX', 'FIRSTY', 'DATACLASS', 'PEAKTABLE', 'DATATYPE', 'DATACLASS']


TEXT_DATA_TABLE = '##XYDATA= (X++(Y..Y))\n'
TEXT_DATA_TABLE_XY = '##XYDATA= (XY..XY)\n'
TEXT_PEAK_AUTO = '$$ === CHEMSPECTRA PEAK TABLE AUTO ===\n'
TEXT_PEAK_EDIT = '$$ === CHEMSPECTRA PEAK TABLE EDIT ===\n'
TEXT_PEAK_TABLE = '##PEAKTABLE= (XY..XY)\n'
TEXT_ORIGINAL_METADATA = '$$ === CHEMSPECTRA ORIGINAL METADATA ===\n'


class BaseComposer:
    def __init__(self, core):
        self.core = core
        self.title = None
        self.meta = None
        self.itgs = []
        self.mpys = []
        self.all_itgs = []
        self.refArea = 1
        self.refShift = 0
        self.prepare_itg_mpy()

    def __header_pk_common(self, category):
        return [
            '##TITLE={}\n'.format(self.title),
            '##JCAMP-DX=5.00\n',
            '##DATA TYPE={}PEAKTABLE\n'.format(self.core.typ),
            '##DATA CLASS=PEAKTABLE\n',
            '##$CSCATEGORY={}\n'.format(category),
            '##$CSTHRESHOLD={}\n'.format(self.core.threshold),
            '##MAXX={}\n'.format(self.core.boundary['x']['max']),
            '##MAXY={}\n'.format(self.core.boundary['y']['max']),
            '##MINX={}\n'.format(self.core.boundary['x']['min']),
            '##MINY={}\n'.format(self.core.boundary['y']['min'])
        ]

    def __create_sample_description(self):
        ref_name = (
            self.core.params['ref_name'] or
            self.core.dic.get('$CSSOLVENTNAME', [''])[0]
        )
        ref_value = (
            self.core.params['ref_value'] or
            self.core.dic.get('$CSSOLVENTVALUE', [''])[0]
        )
        select_x = (
            self.core.params['select_x'] or
            self.core.dic.get('$CSSOLVENTX', [''])[0]
        )
        ref_value = ref_value if is_number(ref_value) else 0
        select_x = select_x if is_number(select_x) else 0
        spl_desc = [
            '##$CSSOLVENTNAME={}\n'.format(ref_name or ''),
            '##$CSSOLVENTVALUE={}\n'.format(ref_value or '0'),
            '##$CSSOLVENTX={}\n'.format(select_x or '0'),
        ]

        detector = self.core.params.get('detector')
        jcamp_idx = self.core.params.get('jcamp_idx')
        
        if detector:
            curves = detector['curves']
            curve_to_update = next((curve for curve in curves if curve.get('curveIdx') == jcamp_idx), None)

            if curve_to_update:
                selected_detector = curve_to_update.get('selectedDetector', {})
            
                if isinstance(selected_detector, dict):
                    name = selected_detector.get('name')
                    
                    if name:
                        spl_desc.append('##$DETECTOR={}\n'.format(name))

        return spl_desc

    def __header_original_metadata(self):
        return [
            '\n',
            TEXT_ORIGINAL_METADATA,
        ]

    def gen_headers_root(self):
        return [
            '##TITLE={}\n'.format(self.title),
            '##JCAMP-DX=5.0\n',
            '##DATA TYPE=LINK\n',
            '##BLOCKS=1\n',  # TBD
            '\n'
        ]

    def generate_original_metadata(self):
        content = self.__header_original_metadata()
        if self.core.dic is None: return content

        for key, value in self.core.dic.items():
            if is_metadata_to_be_ignored(key):
                continue
            str_value = value
            str_key = key
            if isinstance(value, list):
                str_value = ', '.join([str(val) for val in value])

            if str_key.startswith('#'):
                str_key = str_key[1:]

            content.append(
                '###{}= {}\n'.format(str_key.upper(), str_value)
            )
        content.append('\n\n')
        return content
    def gen_ending(self):
        return [
            '##END=\n',
            '\n'
        ]

    def gen_spectrum_orig(self):
        c_spectrum_orig = [
            '##NPOINTS={}\n'.format(self.core.xs.shape[0]),
            TEXT_DATA_TABLE
        ]
        if (self.core.data_format and self.core.data_format == '(XY..XY)'):
            c_spectrum_orig = [
                '##NPOINTS={}\n'.format(self.core.xs.shape[0]),
                TEXT_DATA_TABLE_XY
            ]
        c_spectrum_orig.extend(self.core.datatable)
        return c_spectrum_orig

    def gen_headers_peaktable_auto(self):
        return ['\n', TEXT_PEAK_AUTO] + self.__header_pk_common('AUTO_PEAK')

    def gen_auto_peaktable(self):
        content = [
            '##NPOINTS={}\n'.format(calc_npoints(self.core.auto_peaks)),
            TEXT_PEAK_TABLE
        ]
        if not self.core.auto_peaks:
            return content

        auto_x = self.core.auto_peaks['x']
        auto_y = self.core.auto_peaks['y']
        for i, _ in enumerate(auto_x):
            content.append(
                '{}, {}\n'.format(auto_x[i], auto_y[i])
            )

        return content

    def gen_headers_peaktable_edit(self):
        header = self.__header_pk_common('EDIT_PEAK')
        spl_desc = self.__create_sample_description()

        return ['\n', TEXT_PEAK_EDIT] + header + spl_desc

    def gen_edit_peaktable(self):
        content = [
            '##NPOINTS={}\n'.format(calc_npoints(self.core.edit_peaks)),
            TEXT_PEAK_TABLE
        ]
        if not self.core.edit_peaks:
            return content

        edit_x = self.core.edit_peaks['x']
        edit_y = self.core.edit_peaks['y']
        for i, _ in enumerate(edit_x):
            content.append(
                '{}, {}\n'.format(edit_x[i], edit_y[i])
            )

        return content

    def tf_jcamp(self):
        meta = ''.join(self.meta)
        tf = tempfile.NamedTemporaryFile(suffix='.jdx')
        tf.write(bytes(meta, 'UTF-8'))
        tf.seek(0)
        return tf

    def prepare_itg_mpy(self):
        if not hasattr(self.core, 'params'):
            return
        core_itg = self.core.params.get('integration')
        core_mpy = self.core.params.get('multiplicity')
        if (not core_itg) or (not core_mpy):
            return
        rArea = core_itg.get('refArea') or 1
        rFact = core_itg.get('refFactor') or 1
        self.refArea = float(rFact) / float(rArea)
        self.refShift = core_itg.get('shift') or 0
        # = = = = =
        itg_stack = core_itg.get('stack') or []
        mpy_stack = core_mpy.get('stack') or []

        self.all_itgs = itg_stack
        for itg in itg_stack:
            skip = False
            for mpy in mpy_stack:
                if (itg['xL'] == mpy['xExtent']['xL']) and (itg['xU'] == mpy['xExtent']['xU']):     # pylint: disable=c0301
                    mpy['area'] = itg['area']
                    self.mpys.append(mpy)
                    skip = True
                    break
            if not skip:
                self.itgs.append(itg)

    def gen_integration_info(self):
        if len(self.itgs) > 0:
            table = []
            for itg in self.itgs:
                absoluteArea = 0
                if 'absoluteArea' in itg:
                    absoluteArea = itg['absoluteArea']
                table.extend([
                    '({}, {}, {}, {})\n'.format(
                        itg['xL'] - self.refShift,
                        itg['xU'] - self.refShift,
                        float(itg['area']) * self.refArea,
                        absoluteArea,
                    ),
                ])
            return table
        elif self.core.params['integration'].get('edited'):
            return []
        elif 'stack' in self.core.params['integration']:
            itg_stack = self.core.params['integration']['stack'] or []
            if 'originStack' not in self.core.params['integration']:
                itg_stack = self.core.itg_table

            if len(itg_stack) == 0:
                return []
            if 'stack' in self.core.params['multiplicity']:
                mpy_stack = self.core.params['multiplicity']['stack'] or []

                if len(mpy_stack) > 0:
                    for itg in itg_stack:
                        for mpy in mpy_stack:
                            if (itg['xL'] == mpy['xExtent']['xL']) and (itg['xU'] == mpy['xExtent']['xU']):     # pylint: disable=c0301
                                return []
            return itg_stack
        else:
            return self.core.itg_table

    def gen_mpy_integ_info(self):
        ascii_uppercase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        if len(self.mpys) > 0:
            table = []
            for idx, mpy in enumerate(self.mpys):
                table.extend([
                    '({}, {}, {}, {}, {}, {}, {}, {}{})\n'.format(
                        idx + 1,
                        mpy['xExtent']['xL'] - self.refShift,
                        mpy['xExtent']['xU'] - self.refShift,
                        calc_mpy_center(mpy['peaks'], self.refShift, mpy['mpyType']),   # noqa: E501
                        float(mpy['area']) * self.refArea,
                        idx + 1,
                        mpy['mpyType'],
                        ascii_uppercase[idx],
                        coupling_string(mpy['js']),
                    ),
                ])
            return table
        elif self.core.params['multiplicity'].get('edited'):
            return []
        elif 'stack' in self.core.params['multiplicity']:
            if 'originStack' not in self.core.params['integration']:
                return self.core.mpy_itg_table
            return self.core.params['multiplicity']['stack']
        else:
            return self.core.mpy_itg_table

    def gen_mpy_peaks_info(self):
        if len(self.mpys) > 0:
            table = []
            mpy_stack = self.core.params['multiplicity'].get('stack') or []
            for mk in mpy_stack:
                mk_idx = 0
                for idx, mpy in enumerate(self.mpys):
                    if (mpy['xExtent']['xL'] == mk['xExtent']['xL']) and (mpy['xExtent']['xU'] == mk['xExtent']['xU']):     # pylint: disable=c0301
                        mk_idx = idx + 1
                        break
                for p in mk['peaks']:
                    table.extend([
                        '({}, {}, {})\n'.format(
                            mk_idx,
                            p['x'] - self.refShift,
                            p['y'] - self.refShift,
                        ),
                    ])
            return table
        elif self.core.params['multiplicity'].get('edited'):
            return []
        else:
            return self.core.mpy_pks_table

    def gen_simulation_info(self):
        if len(self.core.simu_peaks) > 0:
            table = []
            for simu_peak in self.core.simu_peaks:
                table.extend(['{}\n'.format(simu_peak)])
            return table

        return []
