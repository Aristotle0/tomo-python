from tomopy.local.utility import read_option, get_gnsrc
from tomopy.local.misfitseism import get_adjs
import sys


if __name__ == '__main__':
    help_string = """
This program is to compute the misfits and adjoint sources.
The type of misfit functionals can be specified. The source can be
specified too, which will be set according to the iteration if not.

--help              show help information
--type=l2           the specified misfit funcitonal, L2 by default
--path=data         the current directory by default
--nsrc=200          No. of source, some source according to iterations by default
--comp=x            specify the component of waveform to calculate the misfits, x by default
--fmin=10           minimium frequency for CWT
--fmax=50           maximium frequency for CWT
--nf=50             number of frequency sampling

type of misfit functionals:

* --type=l2         L2 waveform difference
* --type=wavelet    Continuous Wavelet transform
* --type=cc         cross-correlation (undifined)

for instance:
misfit --help
misfit --type=l2
misfit --type=wavelet --fmin=20 --fmax=50 --nf=50
"""
    option_dict = read_option(sys, help_string, 1, 8)

    working_path = option_dict.setdefault('path', '.')
    nsrc = int(option_dict.setdefault('src', get_gnsrc(working_path)))
    mtype = option_dict.setdefault('type', 'l2')
    comp = option_dict.setdefault('comp', 'x')
    fmin = float(option_dict.setdefault('fmin', None))
    fmax = float(option_dict.setdefault('fmax', None))
    nf = float(option_dict.setdefault('nf', None))

    get_adjs(working_path, nsrc, mtype, comp, fmin, fmax, nf)





