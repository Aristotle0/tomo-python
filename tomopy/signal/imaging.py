# ploting for tomopy
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np

def waveform_plotting(waveform, time, label=None, figname=None):
    """ Plotting waveform data sets

    Parameter
    ---------
    waveform : ndarray
        (n, ) for only one waveform, and n is the number of time samplings
        (2, n) for observed and synthetics
    time : ndarray
        time samplings
    label : sequence
        ('synthetics', 'observed data')
    figname : string
        if given, figure will be saved in the file called figname
    """
    if waveform == 2 and len(label) != 2:
        raise Exception("Two labels need to be given.")
    fig = plt.figure(figsize=(16, 12), dpi=160)
    ax1 = fig.add_subplot(1, 1, 1)

    if waveform.ndim == 1:
        ax1.plot(time, waveform)
    elif waveform.ndim == 2:
        ax1.plot(time, waveform[0, :], label=label[0])
        ax2.plot(time, waveform[1, :], label=label[1])

    ax1.set_xlim((time[0], time[-1]))
    ax1.set_ylim((np.amin(waveform), np.amax(waveform)))
    ax1.set_xlabel('time(seconds)')
    ax1.set_ylabel('amplitude')

    if figname is None:
        plt.show()
    else:
        plt.savefig(figname)
        plt.close('all')

def gather_plotting(gather, time, loc_rec, label=None, figname=None, wp=1.):
    """ Plotting common-shot-gather

    Parameter
    ---------
    gather : common-shot-gather
        (ng, nt) for only one shot
        (2, ng, nt) for observed data and synthetics
    time : ndarray
        time samplings
    loc_rec : ndarray
        location of receivers
    label : sequence
        ('synthetics', 'observed data')
    figname : string
        if given, figure will be saved in the file called figname
    wp : float
        used for weight when plotting
    """
    if gather == 3 and len(label) != 2:
        raise Exception("Two labels need to be given.")
    fig = plt.figure(figsize=(16, 24), dpi=160)
    ax1 = fig.add_subplot(1, 1, 1)

    weight = np.amax(np.abs(gather)) * wp

    if gather.ndim == 2:
        for i in range(gather.shape[0]):
            ax1.plot(time, gather[i, :]/weight+loc[i], 'b')
    elif gather.ndim == 3:
        colors = ['r', 'b']
        legends = []
        for k in range(2):
            for i in range(gather.shape[1]):
                ax1.plot(time, gather[k, i, :]/weight+loc[i], colors[k])
            legends.append(mlines.Line2D([], [], color=colors[k], label=label[k])
        plt.legend(handles=legends)

    ax1.set_xlim((time[0], time[-1]))
    ax1.set_xlabel('time(seconds)')

    if figname is None:
        plt.show()
    else:
        plt.savefig(figname)
        plt.close('all')

def kernel_plotting(kernel, lx, lz, figname=None):
    """ Plotting Frechet Kernel

    Parameter
    ---------
    kernel : ndarray
        (nx, nz) frechet kernel
    lx : ndarray
        axis in x coordinate
    lz : ndarray
        axis in z coordinate
    """
    x, z = np.meshgrid(lx, lz)

    fig = plt.figure(figsize=(16, 12), dpi=160)
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.pcolormesh(x, z, kernel)
    ax1.set_xlim([x.min(), x.max()])
    ax1.set_ylim([y.min(), y.max()])
    ax1.set_xlabel("distance (km)")
    ax1.set_ylabel("depth (km)")
    plt.colorbar()

    if figname is None:
        plt.show()
    else:
        plt.savefig(figname)
        plt.close('all')

def media_plotting(media, lx, lz, vn, figname=None):
    """ Plotting media

    Parameter
    ---------
    kernel : ndarray
        (nx, nz) frechet kernel
    lx : ndarray
        axis in x coordinate
    lz : ndarray
        axis in z coordinate
    vn : string
        'Vp', 'Vs', 'Density'
    """
    x, z = np.meshgrid(lx, lz)
    if vn == 'Vp' or vn == 'Vs':
        str_colorbar = vn + ' ($m/s$)'
    elif vn == 'Density':
        str_colorbar = vn + ' ($kg/m^3$)'

    fig = plt.figure(figsize=(16, 12), dpi=160)
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.pcolormesh(x, z, kernel)
    ax1.set_xlim([x.min(), x.max()])
    ax1.set_ylim([y.min(), y.max()])
    ax1.set_xlabel("distance (km)")
    ax1.set_ylabel("depth (km)")
    cbar = plt.colorbar()
    cbar.set_label(str_colorbar)

    if figname is None:
        plt.show()
    else:
        plt.savefig(figname)
        plt.close('all')

