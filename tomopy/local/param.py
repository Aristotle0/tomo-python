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



