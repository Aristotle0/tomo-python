from tomopy.signal import separation
import numpy as np
from numpy.testing import assert_array_almost_equal_nulp, assert_allclose


def fft_trans(X):
    """ Perform a FFT along time axis
    """
    X_new = np.fft.fft(X, axis=0)
    return X_new

def read_npz(path):
    """ Read the npz file Vx.npz and Vz.npz, then pretreat them.
    """
    xfile = np.load(path + 'Vx.npz')
    zfile = np.load(path + 'Vz.npz')
    terms = ['obs', 'syn']
    ret = {}
    for term in terms:
        x = fft_trans(xfile['Vx_'+term])
        z = fft_trans(zfile['Vz_'+term])
        tmp_co = np.concatenate((x, z))
        ret[term] = tmp_co
    return ret

# class TestMcwbsmf():
#     def test_projection(self):
#         T = np.ones((2, 1))
#         vecs = 2*T
#         start = 0; end = 1
#         T_proj = separation.projection(T, vecs, start, end)
#         assert_allclose(T_proj, T.reshape(T_proj.shape))
# 
#         T = np.ones((2, 1))
#         vecs = np.array([1, 0]).reshape(T.shape)
#         start = 0; end = 1
#         T_proj = separation.projection(T, vecs, start, end)
#         assert T_proj[0] == 1
# 
#     # def test_mc_wbsmf(self):
#     #     Td = read_npz('data/')
#     #     T = Td['syn']
#     #     Ks = 72; Kf = 3; Nf = 400; Nc = 2
#     #     start1 = 0; end1 = 10
#     #     from time import clock
#     #     start = clock()
#     #     x_accept = separation.mc_wbsmf(T, Nc, Ks, Kf, Nf, start1, end1)
#     #     end = clock()
#     #     print("Total cost time: %.4f s" % ((end-start)))
#     #     print(x_accept.shape)
#     #     assert x.accept.shape == T.shape
# 
#     def test_mc_wbsmf_plot(self):
#         s = np.load('tomopy/signal/x.npy')
#         nt, ns = s.shape
#         nt = nt // 2
#         sx = s[:nt, :]
#         sxt = np.fft.ifft(sx, axis=0)
#         xmax = np.amax(sxt)
#         import matplotlib.pyplot as plt
#         for i in range(ns):
#             plt.plot(sxt[:, i]/xmax/5 + i, 'r')
#         plt.show()
#         assert 1