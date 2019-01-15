import tempfile

TEXT_SPECTRUM_ORIG = '$$ === CHEMSPECTRA SPECTRUM ORIG ===\n'
TEXT_SPECTRUM_EDIT = '$$ === CHEMSPECTRA SPECTRUM EDIT ===\n'
TEXT_DATA_TABLE = '##XYDATA= (X++(Y..Y))\n'
TEXT_ASSIGN_AUTO = '$$ === CHEMSPECTRA PEAK ASSIGNMENTS AUTO ===\n'
TEXT_ASSIGN_EDIT = '$$ === CHEMSPECTRA PEAK ASSIGNMENTS EDIT ===\n'
TEXT_PEAK_ASSIGN = '##PEAK ASSIGNMENTS=(XYA)\n'


def calc_npoints(peaks):
    if peaks:
        return len(peaks['x'])
    return 0


def extrac_sp(sp, key):
    query = sp.dic.get(key, '')
    if type(query) is list:
        return query[0]
    return query


def gen_headers_root(sp):
    return [
        '##TITLE={}\n'.format(sp.title),
        '##JCAMP-DX=5.0\n',
        '##DATA TYPE=LINK\n',
        '##BLOCKS=3\n',
        '\n'
    ]


def header_base(sp):
    return [
        '\n',
        TEXT_SPECTRUM_ORIG,
        '##TITLE={}\n'.format(sp.title),
        '##JCAMP-DX=5.00\n',
        '##DATA TYPE={}\n'.format(sp.datatype),
        '##DATA CLASS=XYDATA\n',
        '##ORIGIN={}\n'.format(extrac_sp(sp, 'ORIGIN')),
        '##OWNER={}\n'.format(extrac_sp(sp, 'OWNER')),
    ]


def header_nmr(sp):
    return [
        '##.OBSERVE FREQUENCY={}\n'.format(extrac_sp(sp, '.OBSERVEFREQUENCY')),
        '##.OBSERVE NUCLEUS={}\n'.format(extrac_sp(sp, '.OBSERVENUCLEUS')),
        '##SPECTROMETER/DATA SYSTEM={}\n'.format(extrac_sp(sp, 'SPECTROMETER/DATASYSTEM')),
    ]


def header_params(sp):
    return [
        '##XUNITS={}\n'.format(sp.label['x']),
        '##YUNITS={}\n'.format(sp.label['y']),
        '##XFACTOR={}\n'.format(sp.factor['x']),
        '##YFACTOR={}\n'.format(sp.factor['y']),
        '##FIRSTX={}\n'.format(sp.boundary['x']['max']),
        '##LASTX={}\n'.format(sp.boundary['x']['min']),
        '##MAXX={}\n'.format(sp.boundary['x']['max']),
        '##MAXY={}\n'.format(sp.boundary['y']['max']),
        '##MINX={}\n'.format(sp.boundary['x']['min']),
        '##MINY={}\n'.format(sp.boundary['y']['min'])
    ]


def gen_headers_spectrum_orig(sp):
    if sp.typ == 'INFRARED':
        return header_base(sp) + header_params(sp)
    else:
        return header_base(sp) + header_nmr(sp) + header_params(sp)


# def gen_headers_spectrum_edit(sp):
#     return [
#         '\n',
#         TEXT_SPECTRUM_EDIT,
#         '##TITLE=SPECTRUM_EDIT_{}\n'.format(sp.title),
#         '##JCAMP-DX=5.00\n',
#         '##DATA TYPE={}\n'.format(sp.datatype),
#         '##DATA CLASS=NTUPLES\n',
#         '##MAXX={}\n'.format(sp.boundary['x']['max']),
#         '##MAXY={}\n'.format(sp.boundary['y']['max']),
#         '##MINX={}\n'.format(sp.boundary['x']['min']),
#         '##MINY={}\n'.format(sp.boundary['y']['min']),
#     ]


def header_pk_common(sp):
    return [
        '##TITLE={}\n'.format(sp.title),
        '##JCAMP-DX=5.00\n',
        '##DATA TYPE={} PEAK ASSIGNMENTS\n'.format(sp.typ),
        '##DATA CLASS=ASSIGNMENTS\n',
        '##THRESHOLD={}\n'.format(sp.threshold),
        '##MAXX={}\n'.format(sp.boundary['x']['max']),
        '##MAXY={}\n'.format(sp.boundary['y']['max']),
        '##MINX={}\n'.format(sp.boundary['x']['min']),
        '##MINY={}\n'.format(sp.boundary['y']['min'])
    ]


def gen_headers_peakassignments_auto(sp):
    return ['\n', TEXT_ASSIGN_AUTO] + header_pk_common(sp)


def gen_headers_peakassignments_edit(sp):
    select_x = sp.params['select_x']
    ref_name = sp.params['ref_name']
    ref_value = sp.params['ref_value']

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

    spl_desc = [
        '##SAMPLE DESCRIPTION={}\n'.format(select_x + ref_name + ref_value)
    ]

    return ['\n', TEXT_ASSIGN_EDIT] + header_pk_common(sp) + spl_desc


def gen_ending():
    return [
        '##END=\n',
        '\n'
    ]


def get_peaks_table(peaks):
    content = []
    for idx, p in enumerate(peaks):
        target = '({}, {}, <{}>)'.format(p['x'], p['y'], idx + 1)
        content.append(target)
    return '\n'.join(content)


def gen_spectrum_orig(sp):
    c_spectrum_orig = [
        '##NPOINTS={}\n'.format(sp.y.shape[0]),
        TEXT_DATA_TABLE
    ]
    c_spectrum_orig.extend(sp.datatable)

    return c_spectrum_orig


# def gen_spectrum_edit(sp):
#     c_spectrum_edit = [
#         '##NPOINTS={}\n'.format(sp.y.shape[0]),
#         TEXT_DATA_TABLE
#     ]
#     c_spectrum_edit.extend(sp.datatable)

#     return c_spectrum_edit


def gen_auto_peakassignments(sp):
    c_peakassignments = [
        '##NPOINTS={}\n'.format(calc_npoints(sp.auto_peaks)),
        TEXT_PEAK_ASSIGN
    ]
    if not sp.auto_peaks:
        return c_peakassignments

    auto_x = sp.auto_peaks['x']
    auto_y = sp.auto_peaks['y']
    for i, _ in enumerate(auto_x):
        c_peakassignments.append('({}, {}, <{}>)\n'.format(auto_x[i], auto_y[i], i + 1))

    return c_peakassignments


def gen_edit_peakassignments(sp):
    c_peakassignments = [
        '##NPOINTS={}\n'.format(calc_npoints(sp.edit_peaks)),
        TEXT_PEAK_ASSIGN
    ]
    if not sp.edit_peaks:
        return c_peakassignments

    edit_x = sp.edit_peaks['x']
    edit_y = sp.edit_peaks['y']
    for i, _ in enumerate(edit_x):
        c_peakassignments.append('({}, {}, <{}>)\n'.format(edit_x[i], edit_y[i], i + 1))

    return c_peakassignments


def count_end_idx(lines):
    end_idx = None
    for idx, v in reversed(list(enumerate(lines))):
        if '##END' in v:
            end_idx = idx
            break
    return end_idx


def count_auto_idx(lines):
    auto_idx = -1
    for idx, v in enumerate(lines):
        if TEXT_ASSIGN_AUTO in v:
            auto_idx = idx
            break
    return auto_idx


def build_block_meta(sp):
    meta = []
    meta.extend(gen_headers_root(sp))

    meta.extend(gen_headers_spectrum_orig(sp))
    meta.extend(gen_spectrum_orig(sp))
    meta.extend(gen_ending())

    # meta.extend(gen_headers_spectrum_edit(sp))
    # meta.extend(gen_spectrum_edit(sp))
    # meta.extend(gen_ending())

    meta.extend(gen_headers_peakassignments_auto(sp))
    meta.extend(gen_auto_peakassignments(sp))
    meta.extend(gen_ending())

    meta.extend(gen_headers_peakassignments_edit(sp))
    meta.extend(gen_edit_peakassignments(sp))
    meta.extend(gen_ending())

    meta.extend(gen_ending())
    return meta


def construct_jcamp_temp(meta):
    meta = ''.join(meta)
    tf = tempfile.NamedTemporaryFile(suffix='.jdx')
    with open(tf.name, 'w') as f:
        tf.write(bytes(meta, 'UTF-8'))
        tf.seek(0)
    return tf


def convert_jcamp_temp(sp):
    meta = build_block_meta(sp)
    return construct_jcamp_temp(meta)

