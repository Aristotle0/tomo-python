# Calculate the estimation of the multicomponent wideband spectral matrix,
# including spatial smoothing and frequential smoothing

import numpy as np
from scipy.linalg import svd

def freqcut(X, Nf):
    """ Cut data array along frequency axis
    """
    Nt = X.shape[0]//2
    X_cut = np.concatenate((X[:Nf, :], X[Nt:Nt+Nf, :]), axis=0)
    return X_cut.flatten()

def freqrecover(T, Ns, Nf, Nt):
    """ recover frequency axis to Nt
    """
    T_new = T.reshape(Nf*2, Ns)
    X = np.zeros((Nt*2, Ns), dtype=np.complex)
    X[:Nf, :] = T_new[:Nf, :]
    X[Nt:Nf+Nt, :] = T_new[Nf:, :]
    return X

@profile
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

    N = T.shape[0]//2
    # spatial smoothing
    # indx1 = []
    # for f in range(Nf):
    #     for s in range(ks-1, ks-1+Ns-2*Ks):
    #         indx1.append(f*Ns+s)
    indx1 = [f*Ns+s for f in range(Nf) for s in range(ks-1, ks-1+Ns-2*Ks)]
    indx2 = [x+int(N//2) for x in indx1]
    indx = indx1 + indx2
    mask = np.zeros_like(T, dtype=bool)
    mask[indx] = True
    T_ss = np.where(mask, T, 0.)

    # frequential smoothing
    T_extend = extend_T(T_ss, Ns, Nf, Nc, Kf)
    # indx1 = []
    # for f in range(kf-1, kf-1+Nf):
    #     for s in range(Ns):
    #         indx1.append(f*Ns+s)
    indx1 = [f*Ns+s for f in range(kf-1, kf-1+Nf) for s in range(Ns)]
    indx2 = [x+N for x in indx1]
    indx = indx1 + indx2
    T_fs = T_extend[indx]

    return T_fs.reshape(-1, 1)

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
    #T_extend = extend_T(T, Ns, Nf, Nc, Kf)
    tmp = tuple(select_subarray(T, ks, kf, Ns, Nf, Nc, Ks, Kf)
        for ks in range(1, 2*Ks+2) for kf in range(1, 2*Kf+2))
    C = np.concatenate(tmp, axis=1)
    return C

@profile
def eig(C):
    """ Eigenvalues and Eigenvectors for C^H*C, and sort the eigenvalues
    from the biggest to the smallest.
    """
    tmp = np.dot(C.conj().T, C)
    s, v, d = svd(tmp)
    # vecs = np.zeros_like(C, dtype=np.complex)
    print(s.shape, v.shape)
    vecs = np.dot(C, s/v)
    # for i in range(C.shape[1]):
        # vecs[:, i] = np.dot(C, s[:, i])/v[i]
    return vecs

def projection(T, vecs, start, end):
    """ Return the projection of T onto the given subspace
    """

    T = T.reshape(-1,)
    T_proj_accept = np.zeros_like(T, dtype=np.complex)
    K = vecs.shape[1]
    for k in range(start, end):
        vec = vecs[:, k]
        norm = np.sum(vec**2)
        T_proj_accept += np.dot(T, vec) / norm * vec
    return T_proj_accept

@profile
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
    vecs = eig(C)
    print(vecs.shape)
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
    start1 = 0; end1 = 10
    from time import clock
    start = clock()
    x_accept = mc_wbsmf(T, Nc, Ks, Kf, Nf, start1, end1)
    # T = extend_T(T, Ns, Nf, Nc, Kf)
    # x = select_subarray(T, 1, 1, Ns, Nf, Ks, Kf)
    # x = np.dot(x, x.conj().T)
    end = clock()
    print("Total cost time: %.4f s" % ((end-start)))
    print(x_accept.shape)
    np.save('x.npy', x_accept)

    # from scipy.sparse.linalg import eigs
    # vals, vecs = eigs(x, k=10)
    # import matplotlib.pyplot as plt
    # plt.plot(np.abs(vals))
    # plt.show()

