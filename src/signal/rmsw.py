# Calculate the estimation of the multicomponent wideband spectral matrix,
# including spatial smoothing and frequential smoothing

import numpy as np
from scipy.linalg import svd

def freqcut(X, Nf):
    """ Cut data array along frequency axis
    """
    Nt = X.shape[0]/2
    X_cut = np.concatenate((X[:Nf, :], X[Nt:Nt+Nf, :]), axis=0)
    return X_cut.flatten()

def freqrecover(T, Ns, Nf, Nt):
    """ recover frequency axis to Nt
    """
    T_new = T.reshape(Nf, Ns)
    X = np.zeros((Nt*2, Ns), dtype=np.complex)
    X[:Nf, :] = T[:Nf, :]
    X[Nt:Nf+Nt, :] = T[Nf:, :]
    return X

def get_index(ks, kf, N, Ns, Nf, Ks, Kf):
    """ Get the index of subarray

    N is the size of array T
    """
    indx = []
    for f in range(kf-1, kf-1+Nf):
        for s in range(ks-1, ks-1+Ns-2*Ks):
            indx.append(f*Ns+s)
    indx2 = [x+int(N/2) for x in indx]
    return indx + indx2

def select_subarray(T, ks, kf, Ns, Nf, Ks, Kf):
    """ Select a subarray from the array T

    T: [X_1(f_1), X_2(f_1), ... X_Ns(f_1),
        X_1(f_2), X_2(f_2), ... X_Ns(f_2),
        ...
        X_1(f_Nf), X_2(f_Nf), ... X_Ns(f_Nf),
        ...
        Z_1(f_1), ...
        ...]
    """
    indx = get_index(ks, kf, T.size, Ns, Nf, Ks, Kf)
    T_sub = T[indx] / np.sqrt((2*Ks+1)*(2*Kf+1))
    return T_sub.reshape(-1, 1)

def extend_T(T, Ns, Nf, Nc, Kf):
    """ Extend T and add the ficticius frequential bands

    """
    tmp = T.reshape(Nf*Nc, Ns)
    col_pad = np.zeros((Nf*Nc, Kf))
    tmp = np.hstack((col_pad, tmp, col_pad))
    return tmp.flatten()

def obs_mat(T, Ns, Nf, Nc, Ks, Kf):
    """ Get the observation matrix of size (M x K)
    """
    M = (Ns-2*Ks)*Nf*Nc
    T_extend = extend_T(T, Ns, Nf, Nc, Kf)
    print(T_extend.shape)
    tmp = tuple(select_subarray(T_extend, ks, kf, Ns, Nf, Ks, Kf)
        for ks in range(1, 2*Ks+2) for kf in range(1, 2*Kf+2))
    C = np.concatenate(tmp, axis=1)
    return C

def eig(C):
    """ Eigenvalues and Eigenvectors for C^H*C, and sort the eigenvalues
    from the biggest to the smallest.
    """
    tmp = np.dot(C.conj().T, C)
    s, v, d = svd(tmp)
    vals = np.sqrt(v)
    vecs = np.dot(C, np.dot(vals, v))
    return vals, vecs

def projection(T, vecs, start, end):
    """ Return the projection of T onto the given subspace
    """
    print(T.shape, vecs.shape)
    x = np.linalg.lstsq(T.reshape(-1, 1), vecs)
    proj = np.dot(T, x)
    T_proj_accept = np.sum(proj[:, start:end], axis=1)
    return T_proj_accept

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
    Nt, Ns = X.shape
    Nt = Nt / 2
    T = freqcut(X, Nf)
    C = obs_mat(T, Ns, Nf, Nc, Ks, Kf)
    vals, vecs = eig(C)
    T_accept = projection(T, vecs, start, end)
    X_accept = freqrecover(T_accept, Ns, Nf, Nt)
    if rej == True:
        T_rej = T - T_accept
        X_rej = freqrecover(T_rej, Ns, Nf, Nt)
        return X_accept, X_rej
    else:
        return X_accept

if __name__ == '__main__':
    import npz
    Td = npz.read_npz('../../data/')
    T = Td['syn']
    Ks = 72; Kf = 3; Nf = 400; Nc = 2
    start = 0; end = 10
    from time import clock
    start = clock()
    x_accept = mc_wbsmf(T, Nc, Ks, Kf, Nf, start, end)
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

