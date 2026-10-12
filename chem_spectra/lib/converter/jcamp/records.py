"""The JCAMP-DX records a block declares, block by block.

nmrglue collects every block's records into one list per label, so an index
into that list is a block number only when every block declares the record.
A LINK file whose interferogram declares `##UNITS= CM, VOLTS, ARBITRARY`
beside an absorbance spectrum that declares none is enough to break the
correspondence, and the spectrum is then labelled with the interferogram's
units.

This module rebuilds the blocks from the file itself. It lives apart from
the converters because the file has to be read while it still exists: at the
endpoint the upload is a NamedTemporaryFile that is closed, and so deleted,
before the technique converter is built.
"""

UNIT_RECORDS = ('XUNITS', 'YUNITS', 'UNITS')

# what read_block_records is asked for: the units, plus the datatype that
# says which block each one belongs to.
# The records this app writes itself. They have to be read back from the
# block they were written into: a `$CSINVERTY` on a peak-table block was
# flipping the spectrum's viewport, and a `$CSTRANSMITTANCE` there
# relabelled untouched absorbance as %T.
MANAGED_RECORDS = ('$CSTRANSMITTANCE', '$CSINVERTY')

BLOCK_RECORDS = UNIT_RECORDS + ('DATATYPE',) + MANAGED_RECORDS


def _label_key(label):
    """A JCAMP-DX label as nmrglue keys it: upper case, without spaces,
    dashes, slashes or underscores."""
    return (label.strip().upper().replace(' ', '').replace('-', '')
            .replace('_', '').replace('/', ''))


def read_block_records(path, keys):
    """The `keys` records of each block, in file order, or None.

    Each ##TITLE= opens a block, as JCAMP-DX requires. Only single-line
    values are kept, which is all the unit records need. `utf-8-sig` so a
    byte-order mark cannot hide the first ##TITLE= and shift every block.
    """
    try:
        with open(path, encoding='utf-8-sig', errors='ignore') as handle:
            return _scan(handle, keys)
    except (OSError, TypeError, ValueError):
        return None


# The labels that open a data table. They are not records, so they are not
# collected -- but which block carries the *spectrum* decides which block's
# units describe it, so those are marked.
SPECTRUM_LABELS = frozenset({
    'XYDATA', 'XYPOINTS', 'DATATABLE', 'NTUPLES', 'RADATA',
})

# A peak table is a side table: a file can carry one in its own block, with
# its own units, beside the spectrum. Marking it as the data block is how the
# spectrum came to be labelled from the peak table.
DATA_LABELS = SPECTRUM_LABELS | frozenset({'PEAKTABLE', 'PEAKASSIGNMENTS'})

# Set on the block that opens a spectrum data table. A JCAMP label cannot
# produce this key: `_label_key` removes underscores.
HOLDS_SPECTRUM = '_SPECTRUM'


def _scan(handle, keys):
    """Iterate the handle rather than reading it whole: these files run to
    several megabytes and only the header lines are wanted."""
    blocks = []
    for line in handle:
        if '##' not in line:
            # a data row, or a comment -- and that is nearly every line in
            # the file. Checked before any splitting or stripping.
            continue
        line = line.split('$$', 1)[0].strip()
        # `strip` first: nmrglue strips too, so an indented child block's
        # `##TITLE=` opens a block for it as well. Testing the raw line made
        # those blocks invisible here and nowhere else, and the two readings
        # of the file then disagreed about how many blocks it has.
        if not line.startswith('##') or '=' not in line:
            continue
        label, value = line[2:].split('=', 1)
        key = _label_key(label)
        if key in DATA_LABELS:
            if blocks and key in SPECTRUM_LABELS:
                blocks[-1][HOLDS_SPECTRUM] = True
            continue
        if key == 'TITLE':
            blocks.append({})
        elif blocks and key in keys and value.strip():
            # An empty value is no value, for `##DATA TYPE=` as for the rest.
            # nmrglue warns and drops it, so recording it here was the one
            # thing that could make the two sequences disagree on a file this
            # app had composed itself.
            blocks[-1].setdefault(key, value.strip())
    return blocks
