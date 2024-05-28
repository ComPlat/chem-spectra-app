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

def cal_cyclic_volta_shift_offset(cyclic_data):
    if cyclic_data is None:
        return []
    if isinstance(cyclic_data, dict) ==  False:
        return []
    if 'spectraList' not in cyclic_data:
        return []
    
    spectra_list = cyclic_data['spectraList']
    if isinstance(spectra_list, list) ==  False:
        return []
    
    list_offsets = []
    
    for spectra in spectra_list:
        offset = 0.0
        analysed_data = spectra['list']
        arr_has_ref_value = list(filter(lambda x: x['isRef'] == True,  analysed_data))
        if len(arr_has_ref_value) > 0:
            shift = spectra['shift']
            val = shift['val']
            ref_value = arr_has_ref_value[0]
            e12 = ref_value['e12']
            offset = e12 - val
        list_offsets.append(offset)
    
    return list_offsets

def cal_cyclic_volta_shift_prev_offset_at_index(cyclic_data, index=0):
    if cyclic_data is None:
        return 0.0
    if isinstance(cyclic_data, dict) ==  False:
        return 0.0
    if 'spectraList' not in cyclic_data:
        return 0.0

    spectra_list = cyclic_data['spectraList']
    if isinstance(spectra_list, list) ==  False:
        return 0.0
    if len(spectra_list) <  index:
        return 0.0

    offset = 0.0
    
    if index == len(spectra_list):
      return 0.0
    
    spectra = spectra_list[index]
    hasRefPeak = spectra.get('hasRefPeak', False) == True
    shift = spectra['shift']
    if 'prevValue' in shift:
        offset = shift['prevValue'] if hasRefPeak else -shift['prevValue']

    return offset
