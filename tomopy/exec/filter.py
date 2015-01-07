#!/usr/bin/env python
from tomopy.local.filterseism import seism_filter
from tomopy.local.utility import read_option, get_gnsrc
import sys


if __name__ == '__main__':

    help_string = """
This program perform a filter to waveform.
The source can be specified, and will be set corresponding to
iteration if not. A working path can be specified.

--help              help information
--src=12            specify source number
--path=data         default as the current directory

for instance:
filter
filter --src=200
"""
    option_dict = read_option(sys, help_string, 1, 3)

    working_path = option_dict.setdefault('path', '.')
    gnsrc = int(option_dict.setdefault('src', get_gnsrc(working_path)))


    seism_filter(gnsrc, working_path)
