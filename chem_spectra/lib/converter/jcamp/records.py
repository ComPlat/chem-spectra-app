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
BLOCK_RECORDS = UNIT_RECORDS + ('DATATYPE',)


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
            lines = handle.readlines()
    except (OSError, TypeError, ValueError):
        return None
    blocks = []
    for line in lines:
        line = line.split('$$', 1)[0].strip()
        if not line.startswith('##') or '=' not in line:
            continue
        label, value = line[2:].split('=', 1)
        key = _label_key(label)
        if key == 'TITLE':
            blocks.append({})
        elif blocks and key in keys and value.strip():
            blocks[-1].setdefault(key, value.strip())
    return blocks
