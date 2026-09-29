def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


# Legend entries in the combined image are filenames. Long ones ran across the
# plot and hid the curves (issue #265). 50 keeps the common case intact --
# every fixture name in the repo is well under it -- while bounding the worst.
LEGEND_LABEL_MAX = 50


def shorten_label(label, limit=LEGEND_LABEL_MAX):
    """Trim a legend label, keeping the ends that identify the file.

    The middle is dropped rather than the tail: names in this domain differ by
    their suffix as often as their prefix (`..._13C.jdx` against
    `..._1H.jdx`), so truncating the end alone can make two curves
    indistinguishable in the legend.
    """
    if label is None:
        return label
    label = str(label)
    if len(label) <= limit:
        return label
    if limit <= 3:
        return label[:limit]
    keep = limit - 3
    head = (keep + 1) // 2
    tail = keep - head
    return label[:head] + '...' + (label[-tail:] if tail else '')
