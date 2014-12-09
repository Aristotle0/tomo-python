# Calculate the estimation of the multicomponent wideband spectral matrix,
# including spatial smoothing and frequential smoothing

import numpy as np
from numpy.linalg import norm
from scipy.linalg import svd, eig
from scipy.fftpack import fft, ifft

def fft_trans(X, Nc):
    """ FFT transform
    """
    Nt = X.shape[0]//Nc
    Xf = np.concatenate((fft(X[:Nt, :], axis=0), fft(X[Nt:, :], axis=0)), axis=0)
    return Xf

def reconstruction(Xf, Nc):
    """ IFFT transform
    """
    Nt = Xf.shape[0]//Nc
    X = np.concatenate((ifft(Xf[:Nt, :], axis=0), ifft(Xf[Nt:, :], axis=0)), axis=0)
    return np.real(X)

def freqcut(X, Nc, Nf):
    """ Cut data array along frequency axis
    """
    Nt = X.shape[0]//Nc
    X_cut = np.concatenate((X[:Nf, :], X[Nt:Nt+Nf, :]), axis=0)
    return X_cut.flatten()

def freqrecover(T, Nc, Ns, Nf, Nt):
    """ recover frequency axis to Nt
    """
    T_new = T.reshape(Nf*Nc, Ns)
    X = np.zeros((Nt*Nc, Ns), dtype=np.complex)
    fx = T_new[:Nf, :]
    X[:Nf, :] = fx
    X[Nt-Nf:Nt, :] = fx[::-1, :].conj()
    fz = T_new[Nf:, :]
    X[Nt:Nt+Nf, :] = fz
    X[Nc*Nt-Nf:Nc*Nt, :] = fz[::-1, :].conj()
    return X

#@profile
def select_subarray(T, ks, kf, Ns, Nf, Nc, Ks, Kf):
    """ Select a subarray from the array T

    T: [X_1(f_1), X_2(f_1), ... X_Ns(f_1),
        X_1(f_2), X_2(f_2), ... X_Ns(f_2),
        ...
        X_1(f_Nf), X_2(f_Nf), ... X_Ns(f_Nf),
        ...
        Z_1(f_1), ...
        ...]
    """

    N = Ns * Nf
    # spatial smoothing
    # indx1 = []
    # for f in range(Nf):
    #     for s in range(ks-1, ks-1+Ns-2*Ks):
    #         indx1.append(f*Ns+s)
    indx1 = [f*Ns+s for f in range(Nf) for s in range(ks-1, ks-1+Ns-2*Ks)]
    indx2 = [x+N for x in indx1]
    indx = indx1 + indx2
    mask = np.zeros_like(T, dtype=bool)
    mask[indx] = True
    T_ss = np.where(mask, T, 0.)

    # frequential smoothing
    T_extend = extend_T(T_ss, Ns, Nf, Nc, Kf)
    N = (Nf+2*Kf) * Ns
    # indx1 = []
    # for f in range(kf-1, kf-1+Nf):
    #     for s in range(Ns):
    #         indx1.append(f*Ns+s)
    # indx1 = [f*Ns+s for f in range(kf-1, kf-1+Nf) for s in range(Ns)]
    indx1 = list(range((kf-1)*Ns, (kf-1+Nf)*Ns))
    indx2 = list(range((kf-1)*Ns+N, (kf-1+Nf)*Ns+N))
    # f_range = range(kf-1, kf-1+Nf)
    # s_range = range(Ns)
    # indx1 = list(map((lambda f, s:f*Ns+s), f_range, s_range)) #range(kf-1, kf-1+Nf), range(Ns)))
    # indx2 = list(map((lambda f, s:f*Ns+s+N), f_range, s_range)) #range(kf-1, kf-1+Nf), range(Ns)))
    indx = indx1 + indx2
    T_fs = T_extend[indx]

    return T_fs.reshape(-1, 1)

def extend_T(T, Ns, Nf, Nc, Kf):
    """ Extend T and add the ficticius frequential bands

    """
    tmp = T.reshape(Nf*Nc, Ns)
    row_pad = np.zeros((Kf, Ns))
    tmp1 = np.vstack((row_pad, tmp[:Nf, :], row_pad))
    tmp2 = np.vstack((row_pad, tmp[Nf:2*Nf, :], row_pad))
    tmp = np.vstack((tmp1, tmp2))
    return tmp.flatten()

def obs_mat(T, Ns, Nf, Nc, Ks, Kf):
    """ Get the observation matrix of size (M x K)
    """
    tmp = tuple(select_subarray(T, ks, kf, Ns, Nf, Nc, Ks, Kf)
        for ks in range(1, 2*Ks+2) for kf in range(1, 2*Kf+2))
    C = np.concatenate(tmp, axis=1)
    return C

#@profile
def evd(C):
    """ Eigenvalues and Eigenvectors for C^H*C, and sort the eigenvalues
    from the biggest to the smallest.
    """
    tmp = np.dot(C.conj().T, C)
    # s, v, d = svd(tmp)
    v, s = eig(tmp)
    # vecs = np.zeros_like(C, dtype=np.complex)
    vecs = np.dot(C, s/v)
    # for i in range(C.shape[1]):
        # vecs[:, i] = np.dot(C, s[:, i])/v[i]
    return v, vecs

def projection(T, vecs, start, end):
    """ Return the projection of T onto the given subspace
    """

    T = T.reshape(-1,)
    T_proj_accept = np.zeros_like(T, dtype=np.complex)
    K = vecs.shape[1]
    norm = np.sum(vecs**2, axis=0)
    print(vecs.shape, T.shape, norm.shape)
    pre = np.dot(vecs.T, T)/norm
    T_proj_accept = np.dot(vecs[:, start:end], pre[start:end])
    # for k in range(start, end):
    #     vec = vecs[:, k]
    #     norm = np.sum(vec**2)
    #     T_proj_accept += np.dot(T.T, vec) / norm * vec
    return T_proj_accept

def normalize(X, Nc):
    """ Normalize along time axis
    """
    Nt, Ns = X.shape
    Nt = Nt // Nc
    nm = np.zeros((Nc, Ns))
    for i in range(Nc):
        nt_range = range(i*Nt, (i+1)*Nt)
        nm[i, :] = norm(X[nt_range, :], np.inf, axis=0)
        #X[nt_range, :] /= nm[i, :]
        for k in range(Ns):
            if nm[i, k] == 0:
                nm[i, k] = 1.
            X[nt_range, k] /= nm[i, k]
    return X, nm

def inverse_norm(X, nm):
    """ inverse normalize
    """
    Nt, Ns = X.shape
    Nc = nm.shape[0]
    Nt = Nt // Nc
    for i in range(Nc):
        nt_range = range(i*Nt, (i+1)*Nt)
        for k in range(Ns):
            X[nt_range, k] *= nm[i, k]
    return X


#@profile
def mc_wbsmf(X, Nc, Ks, Kf, Nf, start, end, rej=False):
    """ Multicomponent wideband spectral matrix filtering
    
    Parameters:
    X - data array, which dimensions are (Nt*Nc, Ns)
    Nc - number of components
    Ks - order of spatial smoothing
    Kf - order of frequential smoothing
    Nf - index for cut frequency
    start, end - signal subspace responding to eigenvalues
    rej - calculate rejection part of signal if True

    Results:
    X_new - new data array, which dimensions are (Nt*Nc, Ns-2Ks)
    """
    Xn, nm = normalize(X, Nc)
    Nt, Ns = Xn.shape
    Nt = Nt // Nc
    Xnf = fft_trans(Xn, Nc)
    T = freqcut(Xnf, Nc, Nf)
    C = obs_mat(T, Ns, Nf, Nc, Ks, Kf)
    v, vecs = evd(C)
    np.save('v.npy', v)
    T_accept = projection(T, vecs, start, end)
    T_accept = freqrecover(T_accept, Nc, Ns, Nf, Nt)
    X_accept = reconstruction(T_accept, Nc)
    X_accept = inverse_norm(X_accept, nm)
    if rej == True:
        X_rej = X-X_accept
        return X_accept, X_rej
    else:
        return X_accept

if __name__ == '__main__':
    import npz
    Td = npz.read_npz('../../data/')
    T = Td['obs']
    Ks = 74; Kf = 3; Nf = 400; Nc = 2
    start1 = 0; end1 = (2*Ks+1)*(2*Kf+1)
    from time import clock
    start = clock()
    x_accept = mc_wbsmf(T, Nc, Ks, Kf, Nf, start1, end1)
    # T = extend_T(T, Ns, Nf, Nc, Kf)
    # x = select_subarray(T, 1, 1, Ns, Nf, Ks, Kf)
    # x = np.dot(x, x.conj().T)
    end = clock()
    print("Total cost time: %.4f s" % ((end-start)))
    np.save('x.npy', x_accept)

    # from scipy.sparse.linalg import eigs
    # vals, vecs = eigs(x, k=10)
    # import matplotlib.pyplot as plt
    # plt.plot(np.abs(vals))
    # plt.show()

