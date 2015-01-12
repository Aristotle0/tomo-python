#!/usr/bin/env python
from tomopy.local.param import Fd2dParam, get_numpt, get_sta_coord
from tomopy.local.ioseism import read_seism
from tomopy.local.utility import read_option, get_gnsrc
import numpy as np
import matplotlib.pyplot as plt
import sys

if __name__ == '__main__':

    help_string = """
This program plots a gather of seismograms.
Directory and component need to be specified.
The source can be specified, and will be set corresponding to
iteration if not. 

--help              help information
--src=12            specify source number
--dir=syn_seism     specify the directory that seismograms are saved, syn_seism by default
--path=data         the current directory by default
--comp=x            specify which compontent to display, x by default

for instance:
plot_gather --src=200 --dir=syn_seism --comp=x
"""
    option_dict = read_option(sys, help_string, 2, 5)

    nsrc = int(option_dict.setdefault('src', get_gnsrc(working_path)))
    pnm_seism = option_dict.setdefault('dir', 'seism_syn')
    comp = option_dict.setdefault('comp', 'x')
    working_path = option_dict.setdefault('path', '.')

    if comp == 'x':
        id_comp = 0
    elif comp == 'z':
        id_comp = 1

    para = Fd2dParam(working_path)
    dim1, dim2 = para.dim
    nt = para.nt
    seism_gather, coord_gather = gather_seism(working_path, pnm_seism, dim1, dim2, nsrc, nt)

    smax = np.amax(np.abs(seism_gather))
    fig = plt.figure(figsize=(12, 8))
    ax1 = fig.add_subplot(1, 1, 1)
    for nsta in range(coord_gather.shape[0]):
        ax1.plot(time, seism_gather[nsta]/smax + coord_gather[nsta], 'b')
    ax1.set_xlabel("time (s)")
    ax1.set_ylabel("location of receivers (km)")
    ax1.set_xlim([time.min(), time.max()])
    plt.title('the %s component of seismograms for the source %i' % (comp, nsrc))
    plt.show()