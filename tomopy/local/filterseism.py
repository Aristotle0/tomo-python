from ..signal.filter import butter_filter
from .ioseism import read_seism, write_seism
from .param import Fd2dParam, get_numpt
import numpy as np


def seism_filter(gnsrc, working_path):
    Nd = 2
    para = Fd2dParam(working_path)
    dim1, dim2 = para.dim
    nsrc = para.number_of_moment_source
    nt = para.nt
    stept = para.stept
    pnm_syn = para.pnm_syn
    pnm_obs = para.pnm_obs
    pnm_syn_filter = para.pnm_syn_filter
    pnm_obs_filter = para.pnm_obs_filter
    low = para.low_cut
    high = para.high_cut
    btype = para.type_filter

    if btype == 'bandpass':
        Wn = [low, high]
    elif btype == 'lowpass' or btype == 'highpass':
        Wn = [low, ]
    fs = 1./stept

    pnm_raw = [pnm_syn, pnm_obs]
    pnm_raw = [working_path+'/'+x for x in pnm_raw]
    pnm_filter = [pnm_syn_filter, pnm_obs_filter]
    pnm_filter = [working_path+'/'+x for x in pnm_filter]

    for n_i in range(dim1):
        for n_k in range(dim2):
            num_pt = get_numpt(gnsrc, n_i, n_k, working_path)
            if num_pt == 0:
                continue
            for j_pnm in range(len(pnm_raw)):
                seism, time = read_seism(pnm_raw[j_pnm], gnsrc, n_i, n_k, nt, num_pt)
                seism_filter = np.zeros_like(seism)
                for jd in range(Nd):
                    for nsta in range(num_pt):
                        seism_filter[jd, :, nsta] = butter_filter(seism[jd, :, nsta], Wn, btype, fs)
                write_seism(seism_filter, time, pnm_filter[j_pnm], gnsrc, n_i, n_k, nt, num_pt)

