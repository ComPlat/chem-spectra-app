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

def make_equal_space_1d(base):
    if base.data is None and base.dic['XYPOINTS']:
        base.data = __parse_xy_points(base)

    # base.data type is dict
    if isinstance(base.data, dict):
        return base.data

    # base.data type is array
    data_shape = base.data.shape
    if len(data_shape) == 1:
        return base.data
    elif len(data_shape) == 2:
        [_, ys] = base.data.T
        return ys
        # f = interp1d(xs, ys, kind='cubic')
        # x_max = xs.max()
        # x_min = xs.min()
        # num_pts = __num_pts(base, x_max, x_min)
        # new_xs = np.linspace(x_min, x_max, num=num_pts, endpoint=True)
        # ys = f(new_xs)
        # return ys
    else:
        return base.data
