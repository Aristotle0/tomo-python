#!/usr/bin/env python
from tomopy.local.param import Fd2dParam, get_numpt, get_sta_coord
from tomopy.local.ioseism import read_seism
from tomopy.local.user_exception import IllegalArgumentError
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import sys

if __name__ == '__main__':
    if len(sys.argv) < 2 or len(sys.argv) > 5:
        raise IllegalArgumentError("number of arguments doesn't match.")
    else:
        option = sys.argv[1:]
        option_rdash = [s[2:] for s in option]
        option_dict = {}
        for opn in option_rdash:
            if opn == 'help':
                print("""
This program plots a gather of seismograms.
No. of source, directory and component need to be specified.

--help              help information
--src=12            specify source number
--dir=syn_seism     specify the directory that seismograms are saved
--path=data         default as the current directory
--comp=x            specify which compontent to display

for instance:
plot_gather --src=200 --dir=syn_seism --comp=x
                    """)
                sys.exit()
            else:
                k, v = opn.split('=')
                option_dict[k] = v

    nsrc = int(option_dict['src'])
    pnm_seism = option_dict['dir']
    comp = option_dict['comp']
    working_path = option_dict.setdefault('path', '.')

    if comp == 'x':
        id_comp = 0
    elif comp == 'z':
        id_comp = 1

    pnm_seism = working_path +'/'+ pnm_seism
    para = Fd2dParam(working_path)
    dim1, dim2 = para.dim
    nt = para.nt
    stept = para.stept
    seism_gather = []
    coord_gather = []
    for ni in range(dim1):
        nk = dim2 - 1
        num_pt = get_numpt(nsrc, ni, nk, working_path)
        if (num_pt > 0):
            seism, time = read_seism(pnm_seism, nsrc, ni, nk, nt, num_pt)
            # print('seism_shape: ', seism.shape)
            coordx_sta = get_sta_coord(nsrc, ni, nk, working_path)
            for nsta in range(num_pt):
                seism_gather.append(seism[id_comp, :, nsta])
                coord_gather.append(coordx_sta[nsta])
    seism_gather = np.array(seism_gather)
    coord_gather = np.array(coord_gather)

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