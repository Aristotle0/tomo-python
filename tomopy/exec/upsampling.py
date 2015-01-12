from tomopy.local.param import Fd2dParam
from tomopy.local.utility import read_option
from tomopy.local.inv2media import import_inv, exp_media, myinterp2d, get_coord
import sys


if __name__ == '__main__':

    help_string = """
This program is to upsampling the inversion results.
Because the grid of forward computation is different from
the one of inversion, we need to put coaser grid of inversion
results to finer grid.

--help          show help information
--path=data     current directory by default
    """
    option = read_option(sys, help_string, 2, 2)
    working_path = option_dict.setdefault('path', '.')

    para = Fd2dParam(working_path)
    dim1, dim2 = para.dim

    tmp_path = working_path + '/' + 'temp'
    lam_inv, mu_inv, x_inv, z_inv = import_inv(tmp_path)

    coord_path = working_path + '/' + 'input'
    for n_i in range(dim1):
        for n_k in range(dim2):
            x, z = get_coord(coord_path, n_i, n_k)
            lam = myinterp2d(x_inv, z_inv, lam_inv, x, z)
            mu = myinterp2d(x_inv, z_inv, mu_inv, x, z)
            exp_media(lam, mu, n_i, n_k)
