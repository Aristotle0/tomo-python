from tomopy.local.param import Fd2dParam

def test_Fd2dParam():
    folder = 'data'
    para = Fd2dParam(folder)
    assert para.nt == 4000
    assert para.dim == [16, 4]
    assert para.stept == 0.0045
    assert para.pnm_syn == 'seism_syn'
