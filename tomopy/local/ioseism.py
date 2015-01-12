import numpy as np
# from netCDF4 import Dataset
from scipy.io import netcdf

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
    fnc = netcdf.netcdf_file(filenm, 'r')
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
    fnc = netcdf.netcdf_file(filenm, 'w')#, format='NETCDF3_CLASSIC')
    fnc.createDimension('num_pt', num_pt)
    fnc.createDimension('time', nt)
    for i in range(nd):
        x = fnc.createVariable(Ord[i], 'f', ('time', 'num_pt',))
        x[:, :] = seism[i, :, :]
    t = fnc.createVariable('time', 'f', ('time',))
    fnc.close()


def gather_seism(working_path, pnm_seism, dim1, dim2, nsrc, nt):
    """ gather all seismograms of a source into one array, and location of receivers into another

    Parameters
    ----------
    working_path : string
        path of working directory
    pnm_seism : string
        directory name for gathering seismograms
    dim1, dim2 : int
        number of MPI blocks
    nsrc : int
        No. of source
    nt : float
        number of time samplings

    Returns
    -------
    seism_gather : ndarray
        gathering of seimograms, (nsta, time)
    coord_gather : ndarray
        gathering of receivers' location in x direction. (nsta)
    """
    seism_gather = []
    coord_gather = []
    pnm_seism = working_path +'/'+ pnm_seism
    for ni in range(dim1):
        for nk in range(dim2):
            num_pt = get_numpt(nsrc, ni, nk, working_path)
            if (num_pt > 0):
                seism, time = read_seism(pnm_seism, nsrc, ni, nk, nt, num_pt)
                coordx_sta = get_sta_coord(nsrc, ni, nk, working_path)
                for nsta in range(num_pt):
                    seism_gather.append(seism[id_comp, :, nsta])
                    coord_gather.append(coordx_sta[nsta])
    seism_gather = np.array(seism_gather)
    coord_gather = np.array(coord_gather)
    return seism_gather, coord_gather
