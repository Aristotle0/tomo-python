from tomopy.signal import tf_misfit
from numpy.testing import assert_array_almost_equal_nulp, assert_allclose
class TestTf():
    def test_cwt(self):
        import numpy as np
        from obspy.core import read
        
        st = read()
        tr = st[0]
        npts = tr.stats.npts
        dt = tr.stats.delta
        t = np.linspace(0, dt * npts, npts)
        f_min = 1
        f_max = 50

        print(tr.data.shape)
        cwt = tf_misfit.cwt(tr.data, dt, result='dictionary')
        print('J', cwt['J'])
        icwt = tf_misfit.icwt(cwt['W'], cwt['sj'], cwt['dt'], cwt['dj'])
        
        import pylab as pl
        pl.subplot(211)
        pl.plot(tr.data)
        pl.subplot(212)
        pl.plot(icwt)
        pl.show()
        assert_allclose(tr.data, icwt)
