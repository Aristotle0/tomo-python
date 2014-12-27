import numpy as np
from netCDF4 import Dataset

def read_seism(folder, gnsrc, n_i, n_k, nt, num_pt):
    """ Read seismograms

    Parameters
    ----------
    folder : string
        Folder name where the file we want is
    gnsrc : int
        No. of the seismic source
    n_i, n_k : int
        Id of the MPI block
    nt : int
        Number of time samplings
    num_pt : int
        Number of receivers

    Returns
    -------
    seism : ndarray
        (ndims, nt, num_pt)
    """
    nd = 2
    Ord = ['Vx', 'Vz']
    str1 = (folder, gnsrc, n_i, n_k)
    filenm = '%s/seismo_s%03i_mpi%02i%02i.nc' % str1
    fnc = Dataset(filenm, 'r')
    seism = np.zeros((nd, nt, num_pt))
    for i in range(nd):
        seism[i, :,  :] = fnc.variables[Ord[i]][:, :]
    time = fnc.variables['time'][:]
    fnc.close()
    return seism, time

def write_seism(seism, time, folder, gnsrc, n_i, n_k, nt, num_pt):
    """ Write seismograms

    Parameters
    ----------
    seism : ndarray
        (ndims, nt, num_pt)
    folder : string
        Folder name where the file we want is
    gnsrc : int
        No. of the seismic source
    n_i, n_k : int
        Id of the MPI block
    nt : int
        Number of time samplings
    num_pt : int
        Number of receivers
    """
    nd = 2
    Ord = ['Vx', 'Vz']
    str1 = (folder, gnsrc, n_i, n_k)
    filenm = '%s/seismo_s%03i_mpi%02i%02i.nc' % str1
    fnc = Dataset(filenm, 'w', format='NETCDF3_CLASSIC')
    fnc.createDimension('num_pt', num_pt)
    fnc.createDimension('time', nt)
    for i in range(nd):
        x = fnc.createVariable(Ord[i], 'f', ('time', 'num_pt',))
        x[:, :] = seism[i, :, :]
    t = fnc.createVariable('time', 'f', ('time',))
    fnc.close()


