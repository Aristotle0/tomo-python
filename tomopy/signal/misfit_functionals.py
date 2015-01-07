from .pycwt import *
import numpy as np

__all__ = ['wavemisf', 'waveletmisf']


def wavemisf(obs, syn):
    """ L2 waveform misfit
    """
    misf = sum(1./2*(obs-syn)**2)
    adjs = (syn - obs)[::-1]
    return misf, adjs

def waveletmisf(obs, syn, fmin, fmax, nf=100):
    """ wavelet misfit
    """
    scales = np.logspace(np.log10(fmin), np.log10(fmax), nf)
    mother_wavelet = Morlet(len_signal=len(obs), scales=scales, normalize=True)
    wt_obs = cwt(obs, mother_wavelet).coefs
    wt_syn = cwt(syn, mother_wavelet).coefs

    dphase = np.angle(wt_syn * wt_obs.conj())
    nrm = np.sum(np.abs(wt_syn)**2)
    misf = np.sum(dphase*np.abs(wt_syn)**2)
    misf /= (2*nrm*np.pi**2)

    weight_function = lambda x: x**(-0.5) 

    ker1 = (dphase**2-misf*np.pi*2)*wt_syn
    wavecoefs = Wavelet(ker1, mother_wavelet, weight_function, syn.dtype)
    wave1 = icwt(wavecoefs)

    ker2 = dphase*wt_syn
    wavecoefs = Wavelet(ker2, mother_wavelet, weight_function, syn.dtype)
    wave2 = icwt(wavecoefs)

    adjs = (-np.real(wave1) + np.imag(wave2))
    adjs /= (nrm*np.pi**2)

    return misf, adjs


