import tempfile
from chem_spectra.lib.shared.misc import is_number


def extrac_dic(core, key):
    query = core.dic.get(key, '')
    if type(query) is list:
        return query[0]
    return query


def calc_npoints(peaks):
    if peaks:
        return len(peaks['x'])
    return 0


TEXT_DATA_TABLE = '##XYDATA= (X++(Y..Y))\n'
TEXT_PEAK_AUTO = '$$ === CHEMSPECTRA PEAK TABLE AUTO ===\n'
TEXT_PEAK_EDIT = '$$ === CHEMSPECTRA PEAK TABLE EDIT ===\n'
TEXT_PEAK_TABLE = '##PEAKTABLE= (XY..XY)\n'


class BaseComposer:
    def __init__(self, core):
        self.core = core
        self.title = None
        self.meta = None

    def __header_pk_common(self):
        return [
            '##TITLE={}\n'.format(self.title),
            '##JCAMP-DX=5.00\n',
            '##DATA TYPE={}PEAKTABLE\n'.format(self.core.typ),
            '##DATA CLASS=PEAKTABLE\n',
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

        return spl_desc

    def gen_headers_root(self):
        return [
            '##TITLE={}\n'.format(self.title),
            '##JCAMP-DX=5.0\n',
            '##DATA TYPE=LINK\n',
            '##BLOCKS=1\n',  # TBD
            '\n'
        ]

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
        c_spectrum_orig.extend(self.core.datatable)
        return c_spectrum_orig

    def gen_headers_peaktable_auto(self):
        return ['\n', TEXT_PEAK_AUTO] + self.__header_pk_common()

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
        header = self.__header_pk_common()
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
