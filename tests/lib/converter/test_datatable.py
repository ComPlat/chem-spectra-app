from chem_spectra.lib.converter.datatable import DatatableModel
import numpy as np


def test_datatable_encode_empty():
    data_table_model = DatatableModel()
    calculate_ys = data_table_model.encode(arr_ys=np.array(
        []), y_factor=None, arr_xs=None, is_xypoints=False)
    assert calculate_ys == []


def test_datatable_encode_xypoints():
    data_table_model = DatatableModel()
    arr_xs = np.array([1.49048, 1.48049])
    arr_ys = np.array([5.34724, 4.81545])
    calculate_ys = data_table_model.encode(
        arr_ys=arr_ys, y_factor=None, arr_xs=arr_xs, is_xypoints=True)
    assert calculate_ys == [
        '1.49048, 5.34724\n',
        '1.48049, 4.81545\n']
