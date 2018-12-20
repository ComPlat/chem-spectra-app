import tempfile


TEXT_ASSIGN_AUTO = '$$ === CHEMSPECTRA PEAK ASSIGNMENTS AUTO ===\n'
TEXT_ASSIGN_EDIT = '$$ === CHEMSPECTRA PEAK ASSIGNMENTS EDIT ===\n'
TEXT_PEAK_ASSIGN = '##PEAK ASSIGNMENTS=(XYA)\n'


def calc_npoints(peaks):
    if peaks:
        return len(peaks['x'])
    return 0


def gen_headers_root(sp):
    return [
        '##TITLE={}_ROOT\n'.format(sp.title),
        '##JCAMP-DX=5.0\n',
        '##DATA TYPE=LINK\n',
        '##BLOCKS={}\n'.format(sp.block_count + 2),
        '\n'
    ]


def gen_headers_peakassignments_auto(sp):
    return [
        '\n',
        TEXT_ASSIGN_AUTO,
        '##TITLE={}_PEAK_ASSIGNMENTS\n'.format(sp.title),
        '##JCAMP-DX=5.00\n',
        '##DATA TYPE={} PEAK ASSIGNMENTS\n'.format(sp.typ),
        '##DATA CLASS=ASSIGNMENTS\n',
        '##THRESHOLD={}\n'.format(sp.threshold),
        '##MAXX={}\n'.format(sp.boundary['x']['max']),
        '##MAXY={}\n'.format(sp.boundary['y']['max']),
        '##MINX={}\n'.format(sp.boundary['x']['min']),
        '##MINY={}\n'.format(sp.boundary['y']['min'])
    ]


def gen_headers_peakassignments_edit(sp):
    return [
        '\n',
        TEXT_ASSIGN_EDIT,
        '##TITLE={}_PEAK_ASSIGNMENTS\n'.format(sp.title),
        '##JCAMP-DX=5.00\n',
        '##DATA TYPE={} PEAK ASSIGNMENTS\n'.format(sp.typ),
        '##DATA CLASS=ASSIGNMENTS\n',
        '##THRESHOLD={}\n'.format(sp.threshold),
        '##MAXX={}\n'.format(sp.boundary['x']['max']),
        '##MAXY={}\n'.format(sp.boundary['y']['max']),
        '##MINX={}\n'.format(sp.boundary['x']['min']),
        '##MINY={}\n'.format(sp.boundary['y']['min'])
    ]


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


def gen_origin(path):
    orig = None
    with open(path, 'r', errors = 'ignore') as f:
        orig = f.readlines()
    return orig


def build_one_block_meta(sp, orig):
    meta = []
    meta.extend(gen_headers_root(sp))
    meta.extend(orig)
    meta.extend(gen_headers_peakassignments_auto(sp))
    meta.extend(gen_auto_peakassignments(sp))
    meta.extend(gen_ending())
    meta.extend(gen_headers_peakassignments_edit(sp))
    meta.extend(gen_edit_peakassignments(sp))
    meta.extend(gen_ending())
    meta.extend(gen_ending())
    return meta


def build_block_meta(sp, orig):
    meta = []
    end_idx = count_end_idx(orig)
    auto_idx = count_auto_idx(orig)
    target_idx = auto_idx if auto_idx > 0 else end_idx
    meta.extend(orig[:target_idx])
    meta.extend(gen_headers_peakassignments_auto(sp))
    meta.extend(gen_auto_peakassignments(sp))
    meta.extend(gen_ending())
    meta.extend(gen_headers_peakassignments_edit(sp))
    meta.extend(gen_edit_peakassignments(sp))
    meta.extend(gen_ending())
    meta.extend(gen_ending())
    return meta


def gen_jcamp_meta(sp, orig):
    if sp.block_count == 1:
        return build_one_block_meta(sp, orig)

    return build_block_meta(sp, orig)


def construct_jcamp_temp(meta):
    meta = ''.join(meta)
    tf = tempfile.NamedTemporaryFile(suffix='.jdx')
    with open(tf.name, 'w') as f:
        tf.write(bytes(meta, 'UTF-8'))
        tf.seek(0)
    return tf


def gen_jcamp_temp(sp, orig):
    meta = gen_jcamp_meta(sp, orig)
    return construct_jcamp_temp(meta)


def edit_jcamp_temp(sp, orig):
    meta = gen_jcamp_meta(sp, orig)
    return construct_jcamp_temp(meta)
