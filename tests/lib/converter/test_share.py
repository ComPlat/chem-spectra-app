import numpy as np
from chem_spectra.lib.converter.share import (
  parse_solvent, reduce_pts
)

def test_parse_solvent():
    assert 1==1
    
def test_reduce_pts_when_does_not_have_any_x():
    array_data = []
    data = np.array(array_data)
    reduced_data = reduce_pts(data)
    assert np.array_equal(reduced_data, data)

def test_reduce_pts_when_under_limitation():
    array_data = [[x, x+1] for x in range(0, 3999)]
    data = np.array(array_data)
    reduced_data = reduce_pts(data)
    assert np.array_equal(reduced_data, data)
    
def test_reduce_pts_when_reached_limitation():
    array_data = [[x, x+1] for x in range(0, 4500)]
    data = np.array(array_data)
    expected_array_data = [[x, x+1] for x in range(576, 4500)]
    expected_data = np.array(expected_array_data)
    reduced_data = reduce_pts(data)
    assert np.array_equal(reduced_data, expected_data)
    