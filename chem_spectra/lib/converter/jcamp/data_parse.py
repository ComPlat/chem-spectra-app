import re

from chem_spectra.lib.converter.share import UnconvertibleSpectrum


class UnparsableJcampData(ValueError):
    """nmrglue read the file but produced no usable data array."""


# Data tables whose length the block declares in ##NPOINTS. NTUPLES pages
# carry their own counts and are not checked here.
COUNTED_TABLES = ('XYDATA', 'XYPOINTS', 'PEAKTABLE')
# LDRs that hold exactly one number. A decimal comma in one of these is how a
# non-conforming writer shows itself -- but only as a hint: jcamp-dx's
# doubleformat.jdx writes `##MAXX=100,1` over a conforming table.
SINGLE_NUMBER_LDRS = (
    'FIRSTX', 'LASTX', 'DELTAX', 'XFACTOR', 'YFACTOR', 'FIRSTY',
    'MINX', 'MAXX', 'MINY', 'MAXY',
)
DECIMAL_COMMA = re.compile(r'^\s*[+-]?\d+,\d+([eE][+-]?\d+)?\s*$')


def _point_count(data):
    """Points in a block's data array, or None when it has no single count."""
    if data is None or isinstance(data, dict):
        return None
    if data.ndim == 3 and data.shape[0] == 1:
        return data.shape[1]
    return data.shape[0]


def refuse_if_point_count_disagrees(block, data):
    """Refuse a block whose data does not have the points it declares.

    JCAMP-DX uses the comma only as a delimiter, so a file written with
    decimal commas parses without error into the wrong numbers: each decimal
    splits in two, and the count no longer matches ##NPOINTS. Across 235
    conforming files the two always agree, so a mismatch means the values
    are not the file's, whatever caused it. Refused rather than drawn.
    """
    if block is None or 'NTUPLES' in block.dataclass:
        return
    if not any(block.has(table) for table in COUNTED_TABLES):
        return
    try:
        declared = int(float(str(block.ldr('NPOINTS')).strip()))
    except (TypeError, ValueError):
        return
    parsed = _point_count(data)
    if parsed is None or parsed == declared:
        return
    reason = ('the data table has {} points but the file declares '
              '##NPOINTS={}, so its values cannot be read reliably'
              .format(parsed, declared))
    commas = [
        '##{}={}'.format(ldr, block.ldr(ldr)) for ldr in SINGLE_NUMBER_LDRS
        if block.has(ldr) and DECIMAL_COMMA.match(str(block.ldr(ldr)))
    ]
    if commas:
        reason += ('. Its header writes decimal commas ({}), which JCAMP-DX '
                   'does not allow; export it with decimal points'
                   .format(commas[0]))
    raise UnconvertibleSpectrum(reason)


def make_ms_data_xsys(base):
    """Every mass-spectrum run in the file, in file order.

    An MS file can hold one run per block. The flat read merged every block's
    pages into one `data['real']` list, which happened to produce this same
    sequence; reading per block, the merge has to be done here, deliberately.

    Each run is whatever shape its block declares -- an NTUPLES page, or the
    `(N, 2)` coordinate pairs of a PEAK TABLE -- which is what `MSComposer`
    already consumes.
    """
    runs = []
    for block in _ms_blocks(base):
        data = block.data
        if data is None:
            continue
        refuse_if_point_count_disagrees(block, data)
        if isinstance(data, dict):
            runs.extend(data.get('real') or [])
        elif data.ndim == 3 and data.shape[0] == 1:
            # (1, N, 2) -- one run of coordinate pairs
            runs.append(data[0])
        else:
            runs.append(data)
    return runs or None


def _ms_blocks(base):
    """The blocks holding mass spectra, or just the target block.

    A core that is not a JCAMP file -- `CdfMSConverter` -- has no file to walk,
    so it falls back to the single block it carries.
    """
    jcamp = getattr(base, 'jcamp', None)
    if jcamp is None:
        return [base.target] if getattr(base, 'target', None) else []
    target = getattr(base, 'target', None)
    if target is None:
        return []
    # By technique, not by spelling. `MASS SPECTRUM` and `CONTINUOUS MASS
    # SPECTRUM` both map to MS, and the flat read merged every block's pages
    # whatever each was called; matching the target's exact string dropped
    # the other spelling's runs without a word. A TIC is not a mass spectrum
    # (it maps to LC/MS), so it stays out.
    typ_of = getattr(base, 'typ_of', None)
    typ = getattr(base, 'typ', '')
    if typ_of is not None and typ:
        blocks = [b for b in jcamp if b.datatype and typ_of(b.datatype) == typ]
    else:
        blocks = [b for b in jcamp if b.datatype == target.datatype]
    return blocks or [target]
