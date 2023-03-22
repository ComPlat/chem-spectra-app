import numpy as np
from decimal import Decimal


# def calc_j(mpy, shift):
#     return 0

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
    if len(ps) == 0:
        return 0
    pxs = list(map(lambda p: p['x'], ps))
    pxs.sort()
    cent_idx = int((len(pxs) - 1) / 2)
    if cent_idx < 0:
        return 0
    return pxs[cent_idx] - shift


def calc_mpy_center(ps, shift, typ):
    count = len(ps)
    avgX = sum(list(map(lambda p: p['x'], ps))) / count - shift
    if typ == 't':
        return centerX(ps, shift) if count == 3 else avgX
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


def to_float(target):
    try:
        target = target.replace(',', '.')
        return float(target)
    except:
        return float(target)


def cal_slope(x1, y1, x2, y2):
    if (x1 == x2):
        return 0
    return (y2-y1)/(x2-x1)

def cal_xyIntegration(xs, ys):
    if len(xs) != len(ys):
        return 0
    # REF: https://github.com/mljs/spectra-processing/blob/master/src/xy/xyIntegration.ts
    integration = 0
    for i in range(0, len(xs)-1):
        integration += ((xs[i + 1] - ys[i]) * (xs[i + 1] + ys[i])) / 2
    return integration

def cal_area_multiplicity(xL, xU, data_xs, data_ys):
    y_max, y_min = np.max(data_ys), np.min(data_ys)
    h = y_max - y_min
    
    iL, iU = len(data_xs)-1, 0
    for (idx, point_x) in enumerate(data_xs):
        if (xL <= point_x and point_x <= xU):
            if (iL > idx):
                iL = idx
            if (idx > iU):
                iU = idx

    ks = []
    k = 0
    for y_value in data_ys:
        cy = y_value / y_max
        if (cy > 0.0):
            k += cy
        ks.append(k)
    
    upper_value = Decimal(str(ks[iU]))
    lower_value = Decimal(str(ks[iL]))

    return float(abs(upper_value - lower_value))
