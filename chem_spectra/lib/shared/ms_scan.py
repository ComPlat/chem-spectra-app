"""How strongly a scan shows the target mass: used to pick the default scan."""

import numpy as np

# The windows around the target m/z, in daltons, that count as the target:
# the ion itself and its +1, +23 (Na) and +39 (K) neighbours.
TARGET_OFFSETS = (0, 1, 23, 39)


def scan_ratios(spc, bound_low, bound_high):
    """(ratio, noise_ratio, max_seed) for one scan of (m/z, intensity) rows.

    - max_seed: the strongest peak inside any target window;
    - ratio: max_seed against the strongest peak up to the last window, in %;
    - noise_ratio: the strongest peak beyond the last window, against the same.

    Vectorised; it replaces a per-point Python loop that was repeated in the
    mzML, JCAMP and netCDF MS converters and gives the same three numbers.
    """
    # Kept in the scan's own type: pymzml hands over float32 for some files,
    # and the loop this replaces did its arithmetic in float32 there. A
    # float64 ratio would differ in the 8th digit, and the ratio is compared
    # (`== 100.0`, `<`) to choose the default scan.
    peaks = np.asarray(spc)
    if not np.issubdtype(peaks.dtype, np.floating):
        peaks = peaks.astype(np.float64)
    if peaks.size == 0:
        xs = ys = np.empty(0, dtype=peaks.dtype)
    else:
        peaks = peaks.reshape(-1, peaks.shape[-1])
        xs, ys = peaks[:, 0], peaks[:, 1]

    seed = np.zeros(xs.shape, dtype=bool)
    for offset in TARGET_OFFSETS:
        seed |= (bound_low + offset < xs) & (xs < bound_high + offset)
    last = bound_high + TARGET_OFFSETS[-1]
    base = xs <= last
    beyond = last < xs

    max_base = ys[base].max() if base.any() else 0.1
    max_seed = ys[seed].max() if seed.any() else 0
    max_oorg = ys[beyond].max() if beyond.any() else 0

    ratio = 100 * max_seed / max_base
    noise_ratio = 100 * max_oorg / max_base
    return ratio, noise_ratio, max_seed
