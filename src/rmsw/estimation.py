# Calculate the estimation of the multicomponent wideband spectral matrix,
# including spatial smoothing and frequential smoothing

import numpy as np

def get_index(ks, kf, N, Ns, Nf, Ks, Kf):
    """ Get the index of subarray

    N is the size of array T
    """
    indx = []
    for f in range(kf-1, kf-1+Nf-2*Kf):
        for s in range(ks-1, ks-1+Ns-2*Ks):
            indx.append(f*Ns+s)
    indx2 = [x+N/2 for x in indx]
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
    T_sub_H = T_sub.conj().T
    return np.outer(T_sub, T_sub_H)

def emsm(T, Ns, Nf, Nc, Ks, Kf):
    """ Estimated multicomponent spectral matrix

    Its dimension is (M x M) with M = (Ns - 2Ks)*(Nf - 2Kf)*Nc
    """
    M = (Ns-2*Ks)*(Nf-2*Kf)*Nc
    ret = np.zeros((M, M), dtype=np.complex)
    for kf in range(1, 2*Kf+2):
        for ks in range(1, 2*Ks+2):
            # print('--%02i--%02i' % (ks, kf))
            ret += select_subarray(T, ks, kf, Ns, Nf, Ks, Kf)
    return ret


if __name__ == '__main__':
    import npz
    Td = npz.read_npz('../../data/')
    T = Td['syn']
    Ns = Td['Ns']
    Nf = Td['Nf']
    Ks = 72; Kf = 3; Nc = 2

    from time import clock
    M = (Ns-2*Ks)*(Nf-2*Kf)*Nc
    x = np.zeros((M, M), dtype=np.complex)
    start = clock()
    # x = emsm(T, Ns, Nf, 2, Ks, Kf)
    for j in range(5):
        for i in range(5):
            x += select_subarray(T, i+1, j+1, Ns, Nf, Ks, Kf)
    end = clock()
    print("cost time: %.4f" % ((end-start)))

