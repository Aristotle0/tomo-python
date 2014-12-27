from tomopy.local.ioseism import read_seism, write_seism
from tomopy.local.param import get_numpt
from numpy.testing import assert_array_almost_equal_nulp, assert_allclose

def test_readseism():
    num_pt = get_numpt(200, 12, 3, 'data')
    s, t = read_seism('data/seism_syn', 200, 12, 3, 4000, num_pt)
    write_seism(s, t, 'data/filter', 200, 12, 3, 4000, num_pt)
    s2, t2 = read_seism('data/filter', 200, 12, 3, 4000, num_pt)
    assert_allclose(s, s2)