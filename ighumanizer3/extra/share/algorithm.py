# -*- coding: utf-8 -*-

__author__ = 'mactep'

import os
import numpy as np


def longest_common_substring(src, dst):
    c = np.zeros((len(src), len(dst)), dtype=np.int)
    z = 0
    src_m = None
    dst_m = None
    for i in range(len(src)):
        for j in range(len(dst)):
            if src[i] == dst[j]:
                if i == 0 or j == 0:
                    c[i, j] = 1
                else:
                    c[i, j] = c[i-1, j-1] + 1
                if c[i, j] > z:
                    z = c[i, j]
                if c[i, j] == z:
                    src_m = (i-z+1, i+1)
                    dst_m = (j-z+1, j+1)
            else:
                c[i, j] = 0
    return src_m, dst_m


def get_common_name(directory, suffix):
    filenames = list(filter(lambda f: f.endswith(suffix), os.listdir(directory)))
    if not filenames:
        return ""
    common = filenames[0]
    for i in xrange(len(filenames)):
        while not filenames[i].startswith(common):
            common = common[:-1]
    return common