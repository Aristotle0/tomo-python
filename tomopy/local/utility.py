from netCDF4 import Dataset
import numpy as np

def associate_blocks(path, fname, vname, nsrc, dim1, dim2, nfd=0):
    """ associate all MPI blocks in one array

    Parameter
    ---------
    path : string
        path of the directory which keeps the target files
    fname : string
        target file name
    vname : string
        target variable name
    nsrc : int
        No. of source
    dim1 : int
        number of blocks in x direction
    dim2 : int
        number of blocks in z direction
    nfd : int
        =0 if no FD layers, =3 if num of FD layers is set 3

    Return
    ------
    xzblocks : ndarray
        array after associating
    """
    for n_i in range(dim1):
        for n_k in range(dim2):
            fnm = ("%s/%s_s%03i_mpi%02i%02i" %
                (path, fname, nsrc, n_i, n_k))
            fnc = Dataset(fnm, 'r')
            block = fnc.variables[vname][:, :]
            block = block[nfd:block.shape[0]-nfd, nfd:block.shape[1]-nfd]
            fnc.close()
            if n_k == 0:
                zblocks = block
            else:
                zblocks = np.concatenate((zblocks, block), axis=0)
        if n_i == 0:
            xzblocks = zblocks
        else:
            xzblocks = np.concatenate((xzblocks, zblocks), axis=1)
    return xzblocks