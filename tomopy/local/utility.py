import numpy as np

class IllegalArgumentError(ValueError):
    pass

def read_option(sys, help_string, nmin, nmax):
    if (len(sys.argv) < nmin or len(sys.argv) > nmax):
        raise IllegalArgumentError("number of arguments doesn't match.")
    else:
        option = sys.argv[1:]
        option_rdash = [s[2:] for s in option]
        option_dict = {}
        for opn in option_rdash:
            if opn == 'help':
                print(help_string)
                sys.exit()
            else:
                k, v = opn.split('=')
                option_dict[k] = v
    return option_dict


def get_gnsrc(path):
    """ Get the number of source in current iteration

    Returns
    -------
    ns : int
        No. of source corresponding to random_sources.txt
    """
    status = np.loadtxt(path+'/sgd_status.log')
    if np.ndim(status):
        niter = len(status)
    else:
        niter = 1
    with open(path+'/random_sources.txt') as filein:
        lines = filein.readlines()
        ns = int(lines[niter-1])
    return ns