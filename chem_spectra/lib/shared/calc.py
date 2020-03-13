import numpy as np

def calc_j(mpy, shift):
    return 0

def calc_ks(ys, y_max, h):
    ks = np.zeros_like(ys)
    cs = ys / y_max

    gate = 0.0
    accum = 0
    for idx, val in enumerate(cs):
        if gate < val:
            accum += val * h * 0.01
        ks[idx] = accum
    return ks

def centerX(ps, shift):
    pxs = list(map(lambda p: p['x'], ps))
    pxs.sort()
    cent_idx = int((len(pxs) - 1) / 2)
    if cent_idx < 0:
        return 0
    return pxs[cent_idx] - shift

def calc_mpy_center(ps, shift, typ):
    count = len(ps)
    avgX = sum(list(map(lambda p: p['x'], ps))) / count - shift
    if typ == 's':
        return avgX
    elif typ == 'd':
        return avgX
    elif typ == 't':
        return centerX(ps, shift) if count == 3 else avgX
    elif typ == 'm':
        return avgX
    else:
        return avgX

def get_curve_endpoint(xs, ys, x1, x2):
    iL, iU = 0, 0
    xL, xU = x1, x2
    if x1 > x2:
        xL, xU = x2, x1
    toggle = False
    for idx, val in enumerate(xs):
        if not toggle and (xL <= val and val <= xU):
            toggle = True
            iL = idx
        if toggle and not (xL <= val and val <= xU):
            toggle = False
            iU = idx
    return iL, iU
