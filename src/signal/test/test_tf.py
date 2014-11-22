from tomopy.signal import tf_misfit
from numpy.testing import assert_array_almost_equal_nulp
class TestTf():
    def test_cwt(self):
        import numpy as np
        from obspy.core import read
        from obspy.signal.tf_misfit import cwt
        
        st = read()
        tr = st[0]
        npts = tr.stats.npts
        dt = tr.stats.delta
        t = np.linspace(0, dt * npts, npts)
        f_min = 1
        f_max = 50
        
        cwt_obspy = cwt(tr.data, dt, 8, f_min, f_max)
        cwt_local = tf_misfit.cwt(tr.data, dt, 8, f_min, f_max) 
        assert_array_almost_equal_nulp(cwt_obspy, cwt_local)