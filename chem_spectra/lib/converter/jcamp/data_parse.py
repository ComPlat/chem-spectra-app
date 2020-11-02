import numpy as np
from scipy.interpolate import interp1d


def __num_pts(base, x_max, x_min):
    if base.typ in ['NMR']:
        return 32000
    elif base.typ in ['INFRARED', 'RAMAN']:
        return int(x_max - x_min + 1) * 2
    else:
        return 32000


def __parse_xy_points(base):
    pts = base.dic['XYPOINTS'][0].split('\n')[1:]
    return np.array([[float(p) for p in pt.split(',')]for pt in pts])

def make_ni_data_ys(base, target_idx):
    if base.data is None and base.dic['XYPOINTS']:
        base.data = __parse_xy_points(base)

    # base.data type is dict
    if isinstance(base.data, dict):
        return base.data['real'][target_idx]

    # base.data type is array
    data_shape = base.data.shape
    if len(data_shape) == 1:
        return base.data
    elif len(data_shape) == 2:
        [_, ys] = base.data.T
        return ys
    else:
        return base.data


def make_ms_data_xsys(base):
    if base.data is None and base.dic['XYPOINTS']:
        base.data = [__parse_xy_points(base)]

    # base.data type is dict
    if isinstance(base.data, dict):
        return base.data['real']

    return base.data
