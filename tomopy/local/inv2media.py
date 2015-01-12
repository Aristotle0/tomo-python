import shutil
import numpy as np
# from netCDF4 import netcdf.netcdf_file
from scipy.io import netcdf
from scipy.interpolate import griddata

def import_inv(path):
    """ Import inversion results
    """
    fnc = netcdf.netcdf_file(path + '/' + 'inversion.nc', 'r')
    lam = fnc.variables['lambda'][:, :]
    mu = fnc.variables['mu'][:, :]
    x = fnc.variables['coord_x'][:, :]
    z = fnc.variables['coord_z'][:, :]
    fnc.close()
    return (lam, mu, x, z)

def get_coord(path, n_i, n_k):
    """ Get the coordinate of grid
    """
    fnm = "%s/coord_mpi%02i%02i.nc" % (path, n_i, n_k)
    fnc = netcdf.netcdf_file(fnm, 'r')
    x = fnc.variables['x'][:, :]
    z = fnc.variables['z'][:, :]
    return x, z


def exp_media(dst, src, lam, mu, n_i, n_k):
    """ update media input files after upsamplint and backup old files
    """
    fsrc = "%s/media_mpi%02i%02i.nc" % (src, n_i, n_k)
    fdst = "%s/media_mpi%02i%02i.nc" % (dst, n_i, n_k)
    shutil.copyfile(fsrc, fdst)
    fnc = netcdf.netcdf_file(fsrc, 'a')#, format='NETCDF3_CLASSIC')
    lamid = fnc.variables['lambda']
    muid = fnc.variables['mu']
    lamid[:, :] = np.exp(lam)
    muid[:, :] = np.exp(mu)
    fnc.close()

def myinterp2d(x, y, z, xnew, ynew, method='linear'):
    """ make 2d interpolation according to the axis
    """
    x = np.ravel(x)
    y = np.ravel(y)
    z = np.ravel(z)
    znew = griddata((x, y), z, (xnew, ynew), method=method, fill_value=0.)
    return znew