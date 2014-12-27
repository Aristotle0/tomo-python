#!/usr/bin/env python
from tomopy.local.filterseism import seism_filter
from tomopy.local.status import get_gnsrc
from tomopy.local.user_exception import IllegalArgumentError
import sys


if __name__ == '__main__':
    if len(sys.argv) > 3:
        raise IllegalArgumentError("number of arguments doesn't match.")
    else:
        option = sys.argv[1:]
        option_rdash = [s[2:] for s in option]
        option_dict = {}
        for opn in option_rdash:
            if opn == 'help':
                print("""
This program perform a filter to waveform.
The source can be specified, and will be set corresponding to
iteration if not. A working path can be specified.

--help              help information
--src=12            specify source number
--path=data         default as the current directory

for instance:
filter
filter --src=200
                    """)
                sys.exit()
            else:
                k, v = opn.split('=')
                option_dict[k] = v

    working_path = option_dict.setdefault('path', '.')
    gnsrc = option_dict.setdefault('src', get_gnsrc(working_path))


    seism_filter(gnsrc, working_path)
