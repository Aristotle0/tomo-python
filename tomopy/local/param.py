from netCDF4 import Dataset
import numpy as np

class Fd2dParam():
    def __init__(self, folder):
        fnm = _contpath(folder, 'SeisFD2D.conf')
        self.get_val(fnm, 'ni', 'int')
        self.get_val(fnm, 'nk', 'int')
        self.get_val(fnm, 'nt', 'int')
        self.get_val(fnm, 'stept', 'float')
        self.get_val(fnm, 'dim', 'int2')

        fnm = _contpath(folder, 'SeisSource.conf')
        self.get_val(fnm, 'number_of_moment_source', 'int')

        fnm = _contpath(folder, 'am.conf')
        self.get_val(fnm, 'pnm_syn', 'string')
        self.get_val(fnm, 'pnm_obs', 'string')
        self.get_val(fnm, 'pnm_syn_filter', 'string')
        self.get_val(fnm, 'pnm_obs_filter', 'string')
        self.get_val(fnm, 'type_filter', 'string')
        self.get_val(fnm, 'low_cut', 'float')
        self.get_val(fnm, 'high_cut', 'float')


    def get_val(self, filename, varname, dtype):
        with open(filename) as infile:
            for line in infile:
                if line.startswith(varname):
                    v = line.split('=')[1]
                    if '#' in v:
                        v = v.split("#", 1)[0]

                    setattr(self, varname, _str2(v, dtype))
                    break
            else:
                raise NameError("Can't find " + varname + " in " + filename)

def _str2(v, dtype):
    if dtype == 'int':
        return int(v)
    elif dtype == 'float':
        return float(v)
    elif dtype == 'string':
        return v.strip()
    elif dtype == 'int2':
        return [int(x) for x in v.split()]
    else:
        raise TypeError

def _contpath(folder, fname):
    if folder.endswith('/'):
        return folder+fname
    else:
        return folder+'/'+fname


def get_numpt(gnsrc, n_i, n_k, folder='.'):
    """ Get number of receivers in one block

    Parameter
    ---------
    gnsrc : int
        No. of the seismic source
    n_i, n_k : int
        Id of the MPI block

    Return
    ------
    num_pt : int
        Number of receivers in one block
    """
    str1 = (folder, gnsrc, n_i, n_k)
    filenm = '%s/input/station_s%03i_mpi%02i%02i.nc' % str1
    fnc = Dataset(filenm, 'r')
    num_pt = len(fnc.dimensions['num_pt'])
    fnc.close()
    return num_pt


def get_sta_coord(gnsrc, n_i, n_k, folder='.'):
    """ Get the location of receivers in x direction.

    Parameter
    ---------
    gnsrc : int
        No. of the seismic source
    n_i, n_k : int
        Id of the MPI block

    Return
    ------
    coordx : ndarray
        x coordinates of receivers (km)
    """
    str1 = (folder, gnsrc, n_i, n_k)
    filenm = '%s/input/station_s%03i_mpi%02i%02i.nc' % str1
    fnc = Dataset(filenm, 'r')
    coordx = fnc.variables['coord'][:, 0]/1.e3
    return coordx


def associate_blocks(path, fname, vname, dim1, dim2, nsrc=None, nfd=0):
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
            if nsrc == None:
                fnm = ("%s/%s_mpi%02i%02i.nc" %
                    (path, fname, n_i, n_k))
            else:
                fnm = ("%s/%s_s%03i_mpi%02i%02i.nc" %
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

