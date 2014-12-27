import numpy as np

def tukeywin(N, alpha=0.5):
    """The Tukey window, also called the tapered cosine window.

    Parameter
    ---------
    N: int
        window length
    alpha: float
        alpha = 0 rectangular window; alpha = 1 hanning window

    Result
    ------
    w: ndarray
        the Tukey window
    """
    if alpha <= 0:
        return np.ones(N)
    elif alpha >= 1:
        return np.hanning(N)

    x = np.linspace(0, 1, N)
    w = np.ones(x.shape)

    first_condition = x < (alpha/2)
    w[first_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[first_condition] - alpha/2)))

    # w = 1 for the second condition, which have been implemented

    third_condition = x >= (1 - alpha/2)
    w[third_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[third_condition] - 1 + alpha/2)))

    return w