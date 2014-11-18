# Preconditions for test data Vx.npz and Vz.npz.
# Their structure is (time, nsta)

__all__ = ['read_npz']
import numpy as np

def fft_trans(X):
    """ Perform a FFT along time axis
    """
    X_new = np.fft.fft(X, axis=0)
    return X_new.flatten()

def read_npz(path):
    """ Read the npz file Vx.npz and Vz.npz, then pretreat them.
    """
    xfile = np.load(path + 'Vx.npz')
    zfile = np.load(path + 'Vz.npz')
    terms = ['obs', 'syn']
    ret = {}
    ret['Nf'], ret['Ns'] = xfile['Vx_obs'].shape
    for term in terms:
        x = fft_trans(xfile['Vx_'+term])
        z = fft_trans(zfile['Vz_'+term])
        tmp_co = np.concatenate((x, z))
        ret[term] = tmp_co
    return ret


if __name__ == '__main__':
    Td = read_npz('../../data/')
    print(Td['obs'].shape, Td['obs'].dtype)


