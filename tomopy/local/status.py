import numpy as np

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
