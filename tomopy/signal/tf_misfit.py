"""
Time Frequency Misfit Fuctions based on [Kristekova2006] and
[Kristekova2009].
cwt and icwt code are referred to pycwt

"""
from numpy import (asarray, arange, array, argsort, arctanh, ceil, concatenate, conjugate, cos, diff, exp, intersect1d, isnan, isreal, log, log2, mod, ones, pi, prod, real, round, sort, sqrt, unique, zeros, polyval, nan, ma, floor, interp, loadtxt, savetxt, angle)
from numpy.fft import fft, ifft, fftfreq
from scipy.signal import convolve2d, lfilter


class Morlet:
    """Implements the Morlet wavelet class.
    Note that the input parameters f and f0 are angular frequencies.
    f0 should be more than 0.8 for this function to be correct, its
    default value is f0=6.
    """

    name = 'Morlet'

    def __init__(self, f0=6.0):
        self._set_f0(f0)

    def psi_ft(self, f):
        """Fourier transform of the approximate Morlet wavelet."""
        return (pi ** -.25) * exp(-0.5 * (f - self.f0) ** 2.)

    def psi(self, t):
        """Morlet wavelet as described in Torrence and Compo (1998)."""
        return (pi ** -.25) * exp(1j * self.f0 * t - t ** 2. / 2.)

    def flambda(self):
        """Fourier wavelength as of Torrence and Compo (1998)."""
        return (4 * pi) / (self.f0 + sqrt(2 + self.f0 ** 2))

    def coi(self):
        """e-Folding Time as of Torrence and Compo (1998)."""
        return 1. / sqrt(2.)

    def sup(self):
        """Wavelet support defined by the e-Folding time."""
        return 1. / coi

    def _set_f0(self, f0):
        # Sets the Morlet wave number, the degrees of freedom and the
        # empirically derived factors for the wavelet bases C_{\delta}, \gamma,
        # \delta j_0 (Torrence and Compo, 1998, Table 2)
        self.f0 = f0             # Wave number
        self.dofmin = 2          # Minimum degrees of freedom
        if self.f0 == 6.:
            self.cdelta = 0.776  # Reconstruction factor
            self.gamma = 2.32    # Decorrelation factor for time averaging
            self.deltaj0 = 0.60  # Factor for scale averaging
        else:
            self.cdelta = -1
            self.gamma = -1
            self.deltaj0 = -1
 
def cwt(signal, dt=1., dj=1./12, s0=-1, J=-1, wavelet=Morlet(), result=None):
    """Continuous wavelet transform of the signal at specified scales.
    PARAMETERS
        signal (array like) :
            Input signal array
        dt (float) :
            Sample spacing.
        dj (float, optional) :
            Spacing between discrete scales. Default value is 0.25.
            Smaller values will result in better scale resolution, but
            slower calculation and plot.
        s0 (float, optional) :
            Smallest scale of the wavelet. Default value is 2*dt.
        J (float, optional) :
            Number of scales less one. Scales range from s0 up to
            s0 * 2**(J * dj), which gives a total of (J + 1) scales.
            Default is J = (log2(N*dt/so))/dj.
        wavelet (class, optional) :
            Mother wavelet class. Default is Morlet()
        result (string, optional) :
            If set to 'dictionary' returns the result arrays as itens
            of a dictionary.
    RETURNS
        W (array like) :
            Wavelet transform according to the selected mother wavelet.
            Has (J+1) x N dimensions.
        sj (array like) :
            Vector of scale indices given by sj = s0 * 2**(j * dj),
            j={0, 1, ..., J}.
        freqs (array like) :
            Vector of Fourier frequencies (in 1 / time units) that
            corresponds to the wavelet scales.
        coi (array like) :
            Returns the cone of influence, which is a vector of N
            points containing the maximum Fourier period of useful
            information at that particular time. Periods greater than
            those are subject to edge effects.
        fft (array like) :
            Normalized fast Fourier transform of the input signal.
        fft_freqs (array like):
            Fourier frequencies (in 1/time units) for the calculated
            FFT spectrum.
    EXAMPLE
        mother = wavelet.Morlet(6.)
        wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(var,
            0.25, 0.25, 0.5, 28, mother)
    """
    n0 = len(signal)                              # Original signal length.
    if s0 == -1: s0 = 2 * dt / wavelet.flambda()  # Smallest resolvable scale
    if J == -1: J = int(log2(n0 * dt / s0) / dj)  # Number of scales
    N = 2 ** (int(log2(n0)) + 1)                  # Next higher power of 2.
    signal_ft = fft(signal, N)                    # Signal Fourier transform
    ftfreqs = 2 * pi * fftfreq(N, dt)             # Fourier angular frequencies

    sj = s0 * 2. ** (arange(0, J+1) * dj)         # The scales
    freqs = 1. / (wavelet.flambda() * sj)         # As of Mallat 1999

    # Creates an empty wavlet transform matrix and fills it for every discrete
    # scale using the convolution theorem.
    W = zeros((len(sj), N), 'complex')
    for n, s in enumerate(sj):
        psi_ft_bar = ((s * ftfreqs[1] * N) ** .5 * 
            conjugate(wavelet.psi_ft(s * ftfreqs)))
        W[n, :] = ifft(signal_ft * psi_ft_bar, N)

    # Checks for NaN in transform results and removes them from the scales,
    # frequencies and wavelet transform.
    sel = ~isnan(W).all(axis=1)
    sj = sj[sel]
    freqs = freqs[sel]
    W = W[sel, :]

    # Determines the cone-of-influence. Note that it is returned as a function
    # of time in Fourier periods. Uses triangualr Bartlett window with non-zero
    # end-points.
    coi = (n0 / 2. - abs(arange(0, n0) - (n0 - 1) / 2))
    coi = wavelet.flambda() * wavelet.coi() * dt * coi
    #
    if result == 'dictionary':
        result = dict(
            W = W[:, :n0],
            sj = sj,
            freqs = freqs,
            #period = 1. / freqs,
            coi = coi,
            signal_ft = signal_ft[1:N//2] / N ** 0.5,
            ftfreqs = ftfreqs[1:N//2] / (2. * pi),
            dt = dt,
            dj = dj,
            s0 = s0,
            J = J,
            wavelet = wavelet
        )
        return result
    else:
        return (W[:, :n0], sj, freqs, coi, signal_ft[1:N//2] / N ** 0.5,
                ftfreqs[1:N//2] / (2. * pi))


def icwt(W, sj, dt, dj=0.25, w=Morlet()):
    """Inverse continuous wavelet transform.
    PARAMETERS
        W (array like):
            Wavelet transform, the result of the cwt function.
        sj (array like):
            Vector of scale indices as returned by the cwt function.
        dt (float) :
            Sample spacing.
        dj (float, optional) :
            Spacing between discrete scales as used in the cwt
            function. Default value is 0.25.
        w (class, optional) :
            Mother wavelet class. Default is Morlet()
    RETURNS
        iW (array like) :
            Inverse wavelet transform.
    EXAMPLE
        mother = wavelet.Morlet(6.)
        wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(var,
            0.25, 0.25, 0.5, 28, mother)
        iwave = wavelet.icwt(wave, scales, 0.25, 0.25, mother)
    """
    a, b = W.shape
    c = sj.size
    if a == c:
        sj = (ones([b, 1]) * sj).transpose()
    elif b == c:
        sj = ones([a, 1]) * sj
    else:
        raise Warning('Input array dimensions do not match.')

    # As of Torrence and Compo (1998), eq. (11)
    iW = dj * sqrt(dt) / w.cdelta * w.psi(0) * (real(W) / sj).sum(axis=0)
    return iW

 