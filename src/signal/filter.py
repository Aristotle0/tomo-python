# Filters for signal processing

from scipy.signal import butter, lfilter

def butter_filter(freq_cut, btype, fs, order=5):
    nyq = 0.5 * fs
    
