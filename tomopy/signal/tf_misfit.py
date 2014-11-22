"""
Time Frequency Misfit Fuctions based on [Kristekova2006] and
[Kristekova2009].

"""

import numpy as np

def morletft(s, t, w0, dt):
    """
    Fourier tranformed morlet function.

    Parameters
    ----------
    s : ndarray, 1-D
        Scales
    t : ndarray, 1-D
        time
    w0 : float
        Parameter for the wavelet, tradeoff between time and frequency
    dt : float
        Time step between two samples in st (in seconds)

    Returns
    -------
    wavelet : ndarray
        fourier transformed morlet function    
    """
    p = 0.7511255444649425 # pi**(-0.25)

    s = s.reshape(-1, 1)
    t = s.reshape(1, -1)
    psi = lambda t: p * np.exp(1j * w0 * t) * np.exp(-t ** 2 / 2.)
    psih = psi(-1 * (t - t[-1] / 2.) / s).conj() / np.abs(s) ** .5
    wavelet = np.fft.fft(psih, n=fft, axis=1)
    return wavelet


def cwt(st, dt, w0, fmin, fmax, nf=100, wl='morlet'):
    """
    Continuous Wavelet Transformation in the Freqency Domain

    Parameters
    ----------
    st : ndarray, 1-D
        Time dependent signal
    dt : float
        Time step between two samples in st (in seconds)
    w0 : float
        Parameter for the wavelet, tradeoff between time and frequency
    fmin : float
        Minimum frequency (in Hz)
    fmax : float
        Maximum frequency (in Hz)
    nf : int
        Number of logarithmically spaced frequencies between fmin and fmax
    wl : string
        Wavelet to use, for now only 'morlet' is valid

    Returns
    -------
    cwt : ndarray, complex
        Time frequency representation of st
    """
    st = asarray(st)
    npts = len(st) * 2
    tmax = (npts - 1) * dt
    t = np.linspace(0., tmax, npts)
    f = np.logspace(np.log10(fmin), np.log10(fmax), nf)
    scale = lambda f: w0 / (2 * np.pi * f)
    s = scale(f)

    nfft = util.newpow2(npts)*2
    sf = np.fft.fft(st, n=nfft)

    if wl == 'morlet':
        wft = morletft(s, t, w0, dt)
    else:
        raise ValueError('wavelet type "' + wl + '" not defined.')

    cwt = np.zeros((nf, npts // 2), dtype=np.complex)

    for i in range(nf):
        tminin = int(t[-1] / 2. / dt)
        cwt[i] = np.fft.ifft(sf * wft[i])[tminin:tminin + npts // 2] * dt

    return cwt




