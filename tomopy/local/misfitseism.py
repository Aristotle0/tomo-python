from ..signal.misfit_functionals import wavemisf, waveletmisf
from .param import Fd2dParam, get_numpt
from .ioseism import read_seism
import numpy as np
# from netCDF4 import Dataset
from scipy.io import netcdf


def write_misf(misf, adjs, folder, gnsrc, n_i, n_k, nt, num_pt, nd):
    """ Write misfits

    Parameters
    ----------
    misf : ndarray
        (num_pt,) misfit values
    ajds : ndarray
        (nt, num_pt) adjoint source function
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
    nd : int
        the component id
    """
    str1 = (folder, gnsrc, n_i, n_k)
    filenm = '%s/misf_s%03i_mpi%02i%02i.nc' % str1
    fnc = netcdf.netcdf_file(filenm, 'w')#, format='NETCDF3_CLASSIC')
    fnc.createDimension('geo', 3)
    fnc.createDimension('num_pt', num_pt)
    fnc.createDimension('time', nt)
    misfid = fnc.createVariable('misf', 'f', ('geo', 'num_pt',))
    misfid[nd, :] = misf
    adjsid = fnc.createVariable('adjs', 'f', ('time', 'geo', 'num_pt'))
    adjsid[:, nd, :] = adjs
    fnc.close()



def get_adjs(working_path, nsrc, mtype, comp, fmin, fmax, nf):

    if mtype == 'l2':
        misf_func = wavemisf
    elif mtype == 'wavelet':
        misf_func = lambda x, y : waveletmisf(x, y, fmin, fmax, nf)

    if comp == 'x':
        id_comp = 0
        id_misf = 0
    elif comp == 'z':
        id_comp = 1
        id_misf = 2

    para = Fd2dParam(working_path)
    dim1, dim2 = para.dim
    nt = para.nt
    pnm_seism = para.pnm_syn_filter, para.pnm_obs_filter
    pnm_seism = [working_path + '/' + x for x in pnm_seism]
    pnm_misf = working_path + '/' + 'misf'

    for n_i in range(dim1):
        for n_k in range(dim2):
            num_pt = get_numpt(nsrc, n_i, n_k, working_path)
            if num_pt == 0:
                continue
            seism_syn, time = read_seism(pnm_seism[0], nsrc, n_i, n_k, nt, num_pt)
            seism_obs, time = read_seism(pnm_seism[1], nsrc, n_i, n_k, nt, num_pt)
            misf = np.zeros(seism_syn.shape[2])
            adjs = np.zeros_like(seism_syn)
            for nsta in range(num_pt):
                misf[nsta], adjs[id_comp, :, nsta] = misf_func(seism_obs[id_comp, :, nsta],
                    seism_syn[id_comp, :, nsta])
            write_misf(misf, adjs[id_comp, :, :], pnm_misf, nsrc, n_i, n_k, nt, num_pt, id_misf)


