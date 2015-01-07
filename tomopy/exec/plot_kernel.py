#!/usr/bin/env python
from tomopy.local.param import Fd2dParam
from tomopy.local.param import associate_blocks
from tomopy.local.utility import read_option, get_gnsrc
import numpy as np
import matplotlib.pyplot as plt
import sys



if __name__ == '__main__':

    help_string = """
This program is to plot the kernel pattern.
The source can be specified, and will be set corresponding to
iteration if not. A working path can be specified.

--help              help information
--src=12            specify source number, default corresponding to iteration
--path=data         default as the current directory
--display=ka        choose ka or kb to display, default ka
--dx=0.025          specify the length of grid in x coordinate
--dz=0.025          default = dx
--cut=0.5           specify display value scope, full value for 1 and half for 0.5

for instance:
plot_kernel
plot_kernel --src=200
plot_kernel --src=200 --path=data --display=ka --dx=0.025 --cut=0.1
"""
    option_dict = read_option(sys, help_string, 1, 6)

    working_path = option_dict.setdefault('path', '.')
    nsrc = int(option_dict.setdefault('src', get_gnsrc(working_path)))
    ker_comp = option_dict.setdefault('display', 'ka')
    dx = float(option_dict.setdefault('dx', 1))
    dz = float(option_dict.setdefault('dz', dx))
    cut = float(option_dict.setdefault('cut', 1))

    para = Fd2dParam(working_path)
    dim1, dim2 = para.dim

    path = working_path +'/' + 'kernel'
    ker = {}
    ker[ker_comp] = associate_blocks(path, 'kernel', ker_comp, dim1, dim2, nsrc=nsrc)

    nz, nx = ker[ker_comp].shape
    gridx = np.arange(0., nx*dx, dx)
    gridz = np.arange(-nz*dz+dz, dz, dz)
    x, z = np.meshgrid(gridx, gridz)

    size_fig = 16
    fig = plt.figure(figsize=(size_fig, size_fig*nz/nx))
    ax1 = fig.add_subplot(1, 1, 1)
    kv = ker[ker_comp]
    pc_min = kv.min()*cut
    pc_max = kv.max()*cut
    img = ax1.pcolormesh(x, z, kv, vmin=pc_min, vmax=pc_max)
    ax1.set_xlim([gridx.min(), gridx.max()])
    ax1.set_ylim([gridz.min(), gridz.max()])
    cb = plt.colorbar(img, format='%.1e', ticks=[pc_min, (pc_min+pc_max)/2, pc_max])
    plt.show()
