import tempfile


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
            '##THRESHOLD={}\n'.format(self.core.threshold),
            '##MAXX={}\n'.format(self.core.boundary['x']['max']),
            '##MAXY={}\n'.format(self.core.boundary['y']['max']),
            '##MINX={}\n'.format(self.core.boundary['x']['min']),
            '##MINY={}\n'.format(self.core.boundary['y']['min'])
        ]


    def __inherit_sample_description(self):
        try:
            target = self.core.dic['SAMPLEDESCRIPTION'][-1]
            if target and (target != ''):
                return target
        except:
            pass

        return ''


    def __create_sample_description(self):
        select_x = self.core.params['select_x']
        ref_name = self.core.params['ref_name']
        ref_value = self.core.params['ref_value']

        if select_x:
            select_x = 'SELECTX={};'.format(select_x)
        else:
            select_x = ''
        if ref_name:
            ref_name = 'SOLNAME={};'.format(ref_name)
        else:
            ref_name = ''
        if ref_value:
            ref_value = 'SOLVAL={};'.format(ref_value)
        else:
            ref_value = ''

        description = select_x + ref_name + ref_value

        if description == '':
            description = self.__inherit_sample_description()

        spl_desc = [
            '##SAMPLE DESCRIPTION={}\n'.format(description)
        ]

        return spl_desc


    def gen_headers_root(self):
        return [
            '##TITLE={}\n'.format(self.title),
            '##JCAMP-DX=5.0\n',
            '##DATA TYPE=LINK\n',
            '##BLOCKS=1\n', # TBD
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
        with open(tf.name, 'w') as f:
            tf.write(bytes(meta, 'UTF-8'))
            tf.seek(0)
        return tf
