from chem_spectra.lib.shared.calc import *
import pytest

# def test_calc_j():
#     #TODO: implement later
#     pass

def test_calc_ks():
    ys = np.array([1.0, 2.0, 3.0, 4.0])
    y_max = 3.0
    h = 3.0
    ks = calc_ks(ys, y_max, h)
    assert np.alltrue(ks == [0.01, 0.03, 0.06, 0.1])

def test_centerX_no_values():
    ps = []
    center_x = centerX(ps=ps, shift=0)
    assert center_x == 0.0

def test_centerX_odd_number_values():
    ps = [{'x':1.0}, {'x':2.0}, {'x':3.0}]
    center_x = centerX(ps=ps, shift=0)
    assert center_x == 2.0

def test_centerX_even_number_values():
    ps = [{'x':1.0}, {'x':2.0}, {'x':3.0}, {'x':4.0}]
    center_x = centerX(ps=ps, shift=0)
    assert center_x == 2.0

def test_centerX_with_shift():
    ps = [{'x':1.0}, {'x':2.0}, {'x':3.0}, {'x':4.0}]
    center_x = centerX(ps=ps, shift=0.5)
    assert center_x == 1.5

def test_calc_mpy_center_not_t_typ():
    ps = [{'x':1.0}, {'x':2.0}, {'x':3.0}, {'x':4.0}]
    mpy_typ = ['s', 'd', 'm']
    for typ in mpy_typ:    
        center_x = calc_mpy_center(ps=ps, shift=0, typ=typ)
        assert center_x == 2.5

def test_calc_mpy_center_t_typ():
    ps = [{'x':1.0}, {'x':2.0}, {'x':3.0}, {'x':4.0}]
    center_x = calc_mpy_center(ps=ps, shift=0, typ='t')
    assert center_x == 2.5

def test_calc_mpy_center_t_typ_three_values():
    ps = [{'x':1.0}, {'x':2.0}, {'x':4.0}]
    center_x = calc_mpy_center(ps=ps, shift=0, typ='t')
    assert center_x == 2.0

def test_get_curve_endpoint():
    #TODO: implement later
    pass

def test_to_float_with_int_value():
    float_value = to_float(10)
    assert type(float_value) is float
    assert float_value == 10.0

def test_to_float_with_float_value():
    float_value = to_float(10.0)
    assert type(float_value) is float
    assert float_value == 10.0

def test_to_float_with_string_value_with_dot():
    float_value = to_float('10.2')
    assert type(float_value) is float
    assert float_value == 10.2

def test_to_float_with_string_value_with_comma():
    float_value = to_float('10,2')
    assert type(float_value) is float
    assert float_value == 10.2

def test_cal_slope_no_slope():
    slope_1 = cal_slope(x1=1.0, x2=1.0, y1=1.0, y2=2.0)
    assert slope_1 == 0
    slope_2 = cal_slope(x1=1.0, x2=2.0, y1=1.0, y2=1.0)
    assert slope_2 == 0

def test_cal_slope():
    slope = cal_slope(x1=1.0, x2=2.0, y1=1.0, y2=2.0)
    assert slope == 1.0

def test_cal_xyIntegration_wrong_value():
    xs = [1.0, 2.0, 3.0]
    ys = [1.0, 2.0, 3.0, 4.0]
    integration_value = cal_xyIntegration(xs=xs, ys=ys)
    assert integration_value == 0

def test_cal_xyIntegration():
    xs = [1.0, 2.0, 3.0]
    ys = [1.0, 2.0, 3.0]
    integration_value = cal_xyIntegration(xs=xs, ys=ys)
    assert integration_value == 4.0

def test_cal_area_multiplicity():
    xL = 1.0
    xU = 2.0
    data_xs = [1.0, 2.0]
    data_ys = [2.0, 4.0]
    area = cal_area_multiplicity(xL, xU, data_xs=data_xs, data_ys=data_ys)
    assert area == 1.0
    
@pytest.fixture
def cyclic_data():
    cyclic_data = {'spectraList':[{'list':[{'min':{'x':-1.48904,'y':-1.10686e-05},'max':{'x':1.80895,'y':9.51171e-06},'isRef':True,'e12':0.15995500000000007,'pecker':{'x':-0.129871,'y':7.78418e-07}}],'selectedIdx':0,'isWorkMaxPeak':True,'jcampIdx':0,'shift':{'ref':None,'val':0,'prevValue': 2.5,}, 'hasRefPeak':True},{'list':[{'min':{'x':-1.48904,'y':-3.3747399999999995e-05},'max':{'x':0.929483,'y':0.00023741},'isRef':True,'e12':-0.27977849999999993}],'selectedIdx':0,'isWorkMaxPeak':True,'jcampIdx':1,'shift':{'ref':None,'val':5}},{'list':[{'min':{'x':0.45977,'y':-0.000226347},'max':{'x':1.00943,'y':0.000371349},'isRef':False,'e12':0.7346}],'selectedIdx':0,'isWorkMaxPeak':True,'jcampIdx':2,'shift':{'ref':None,'val':0}}]}
    return cyclic_data
  
@pytest.fixture
def cyclic_data_no_ref_value():
    cyclic_data = {'spectraList':[{'list':[{'min':{'x':-1.48904,'y':-1.10686e-05},'max':{'x':1.80895,'y':9.51171e-06},'isRef':True,'e12':0.15995500000000007,'pecker':{'x':-0.129871,'y':7.78418e-07}}],'selectedIdx':0,'isWorkMaxPeak':True,'jcampIdx':0,'shift':{'ref':None,'val':0,'prevValue': 2.5,}, 'hasRefPeak':False},{'list':[{'min':{'x':-1.48904,'y':-3.3747399999999995e-05},'max':{'x':0.929483,'y':0.00023741},'isRef':True,'e12':-0.27977849999999993}],'selectedIdx':0,'isWorkMaxPeak':True,'jcampIdx':1,'shift':{'ref':None,'val':5}},{'list':[{'min':{'x':0.45977,'y':-0.000226347},'max':{'x':1.00943,'y':0.000371349},'isRef':False,'e12':0.7346}],'selectedIdx':0,'isWorkMaxPeak':True,'jcampIdx':2,'shift':{'ref':None,'val':0}}]}
    return cyclic_data

def test_cal_cyclic_volta_shift_prev_offset_at_index(cyclic_data):
    expected_offset = 2.5
    offset = cal_cyclic_volta_shift_prev_offset_at_index(cyclic_data, 0)
    assert offset == expected_offset
    
def test_cal_cyclic_volta_shift_prev_offset_at_index_no_ref_value(cyclic_data_no_ref_value):
    expected_offset = -2.5
    offset = cal_cyclic_volta_shift_prev_offset_at_index(cyclic_data_no_ref_value, 0)
    assert offset == expected_offset
