class UnparsableJcampData(ValueError):
    """nmrglue read the file but produced no usable data array."""


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
    blocks = [b for b in jcamp if b.datatype == target.datatype]
    return blocks or [target]
