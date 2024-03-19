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
    pts = []
    if 'XYPOINTS' in base.dic:
        pts = base.dic['XYPOINTS'][0].split('\n')[1:]
    elif base.data_format and base.data_format == '(XY..XY)':
        pts = base.dic['XYDATA_OLD'][0].split('\n')[1:]
    return np.array([[float(p) for p in pt.split(',')]for pt in pts])


def make_ni_data_ys(base, target_idx):
    if base.data is None and base.dic['XYPOINTS']:
        base.data = __parse_xy_points(base)
    elif base.data_format and base.data_format == '(XY..XY)':
        base.data = __parse_xy_points(base)

    # base.data type is dict
    if isinstance(base.data, dict):
        return base.data['real'][target_idx]

    if isinstance(base.data, np.ndarray) == False:
        base.data = np.asarray(base.data)

    # base.data type is array
    data_shape = base.data.shape
    if len(data_shape) == 1:
        return base.data
    elif len(data_shape) == 2:
        try:
            [_, ys] = base.data.T
        except:
            [_, ys] = base.data
        return ys
    else:
        return base.data


def make_ni_data_xs(base):
    if base.data_format and (base.data_format == '(XY..XY)'):
        data = __parse_xy_points(base)
        # base.data type is array
        data_shape = data.shape
        if len(data_shape) == 2:
            [xs, _] = data.T
            return xs

    return None


def make_ms_data_xsys(base):
    if base.data is None and base.dic['XYPOINTS']:
        base.data = [__parse_xy_points(base)]

    # base.data type is dict
    if isinstance(base.data, dict):
        return base.data['real']

    return base.data
