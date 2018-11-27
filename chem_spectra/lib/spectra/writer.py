import tempfile

def gen_headers_root(sp):
    return [
        '##TITLE={}_ROOT\n'.format(sp.title),
        '##JCAMP-DX=5.0\n',
        '##DATA TYPE=LINK\n',
        '##BLOCKS={}\n'.format(sp.block_count + 2),
        '\n'
    ]


def gen_headers_peakassignments_all(sp):
    return [
        '\n',
        '$$ === CHEMSPECTRA PEAK ASSIGNMENTS ALL ===\n',
        '##TITLE={}_PEAK_ASSIGNMENTS\n'.format(sp.title),
        '##JCAMP-DX=5.00\n',
        '##DATA TYPE=PEAK ASSIGNMENTS\n',
        '##DATA CLASS=ASSIGNMENTS\n',
        '##NPOINTS={}\n'.format(sp.peak_idxs.shape[0]),
        '##THRESHOLD={}\n'.format(sp.threshold),
        '##MAXX={}\n'.format(sp.boundary['x']['max']),
        '##MAXY={}\n'.format(sp.boundary['y']['max']),
        '##MINX={}\n'.format(sp.boundary['x']['min']),
        '##MINY={}\n'.format(sp.boundary['y']['min']),
        '##PEAK ASSIGNMENTS=(XYA)\n'
    ]


def gen_headers_peakassignments_edit(sp):
    return [
        '\n',
        '$$ === CHEMSPECTRA PEAK ASSIGNMENTS EDIT ===\n',
        '##TITLE={}_PEAK_ASSIGNMENTS\n'.format(sp.title),
        '##JCAMP-DX=5.00\n',
        '##DATA TYPE=PEAK ASSIGNMENTS\n',
        '##DATA CLASS=ASSIGNMENTS\n',
        '##THRESHOLD={}\n'.format(sp.threshold),
        '##MAXX={}\n'.format(sp.boundary['x']['max']),
        '##MAXY={}\n'.format(sp.boundary['y']['max']),
        '##MINX={}\n'.format(sp.boundary['x']['min']),
        '##MINY={}\n'.format(sp.boundary['y']['min']),
        '$$ --- CHEMSPECTRA EDIT CONTENT ---\n'
    ]


def gen_ending():
    return [
        '##END=\n',
        '\n'
    ]


def gen_content_peakassignments(sp):
    c_peakassignments = []
    for i, idx in enumerate(sp.peak_idxs):
        c_peakassignments.append('({}, {}, <{}>)\n'.format(sp.x[idx], sp.y[idx], i+1))

    c_peakassignments.extend(gen_ending())
    return c_peakassignments


def count_end_idx(lines):
    end_idx = None
    for idx, v in reversed(list(enumerate(lines))):
        if '##END' in v:
            end_idx = idx
            break
    return end_idx


def gen_origin(path):
    orig = None
    with open(path, 'r', errors = 'ignore') as f:
        orig = f.readlines()
    return orig


def gen_jcamp_meta(sp, orig):
    meta = []
    if sp.block_count == 1:
        meta.extend(gen_headers_root(sp))
        meta.extend(orig)
        meta.extend(gen_headers_peakassignments_all(sp))
        meta.extend(gen_content_peakassignments(sp))
        meta.extend(gen_headers_peakassignments_edit(sp))
        meta.extend(gen_ending())
        meta.extend(gen_ending())
    else:
        end_idx = count_end_idx(orig)
        meta.extend(orig[:end_idx])
        meta.extend(gen_headers_peakassignments_all(sp))
        meta.extend(gen_content_peakassignments(sp))
        meta.extend(gen_headers_peakassignments_edit(sp))
        meta.extend(gen_ending())
        meta.extend(gen_ending())
    return meta


def gen_jcamp_temp(sp, orig):
    meta = gen_jcamp_meta(sp, orig)
    meta = ''.join(meta)
    tf = tempfile.NamedTemporaryFile(suffix='.jdx')
    with open(tf.name, 'w') as f:
        tf.write(bytes(meta, 'UTF-8'))
        tf.seek(0)
    return tf
