#!/usr/bin/env python
# Author: Duy Tin Truong (duytin.truong@unitn.it)
#		at CIBIO, University of Trento, Italy

import os
import sys


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep) + sys.path:
            path = path.strip('"')
            exe_file = os.path.join(path, program)

            # print(exe_file)
            if is_exe(exe_file):
                return exe_file

    return None


def is_exe(program):
    return which(program) != None
