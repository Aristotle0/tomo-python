#!/usr/bin/env python
from tomopy.local.param import Fd2dParam
from tomopy.local.param import associate_blocks
from tomopy.local.utility import read_option, get_gnsrc
import numpy as np
import matplotlib.pyplot as plt
import sys

if __name__ == "__main__":

    help_string = """
This program is to plot the media information, such as Vp, Vs and density.
The type of media information needs to be specified.

--help              show help information
--display=vp        specify the media type to display
--path=data         current directory by default 

for instance:
plot_media --help
plot_media --display=vp
"""
    option_dict = read_option(sys, help_string, 2, 3)

    working_path = option_dict.setdefault('path', '.')
    media_comp = option_dict.setdefault('display', 'vp')

    para = Fd2dParam(working_path)
    dim1, dim2 = para.dim

    path = working_path +'/' + 'input'
    media = {}
    keys = ['lambda', 'mu', 'rho']
    for key in keys:
        media[key] = associate_blocks(path, 'media', key, dim1, dim2, nfd=3)

    x = associate_blocks(path, 'coord', 'x', dim1, dim2, nfd=3)
    z = associate_blocks(path, 'coord', 'z', dim1, dim2, nfd=3)
    nz, nx = x.shape

    if media_comp == 'vp':
        media[media_comp] = np.sqrt((media['lambda']+2.*media['mu'])/media['rho'])
        clabel = "Vp (m/s)"
    elif media_comp == 'vs':
        media[media_comp] = np.sqrt(media['mu']/media['rho'])
        clabel = "Vs (m/s)"
    elif media_comp == 'rho':
        clabel = "Density (kg/$m^3)"


    size_fig = 16
    fig = plt.figure(figsize=(size_fig, size_fig*nz/nx))
    ax1 = fig.add_subplot(1, 1, 1)
    mv = media[media_comp]
    pc_min = mv.min()
    pc_max = mv.max()
    img = ax1.pcolormesh(x, z, mv)
    ax1.set_xlim([x.min(), x.max()])
    ax1.set_ylim([z.min(), z.max()])
    ax1.set_xlabel("distance (km)")
    ax1.set_xlabel("depth (km)")
    cb = plt.colorbar(img, format='%.1e', ticks=[pc_min, (pc_min+pc_max)/2, pc_max])
    cb.set_label(clabel)
    plt.show()




