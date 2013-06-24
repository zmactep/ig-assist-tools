# -*- coding: utf-8 -*-

__author__ = 'mactep'

import os
from ighumanizer3.extra.share.fasta_tools import FastaDict, read_fasta
from ighumanizer3.extra.share.alignment import *

CHAIN_MARKER_VL = "VL"
CHAIN_MARKER_VH = "VH"


def load_directory(directory):
    filenames = list(filter(lambda f: f.endswith("-amino.fa"), os.listdir(directory)))
    heavy_chains = FastaDict()
    light_chains = FastaDict()
    for filename in filenames:
        name = os.path.join(directory, filename)
        fd = read_fasta(name)
        for seq in fd.keys():
            if seq.endswith(CHAIN_MARKER_VL):
                light_chains.set(seq, fd.get(seq))
            elif seq.endswith(CHAIN_MARKER_VH):
                heavy_chains.set(seq, fd.get(seq))

    light_chains = multiple_alignment(light_chains, TYPE_UNI_FAST)
    heavy_chains = multiple_alignment(heavy_chains, TYPE_UNI_FAST)

    return light_chains, heavy_chains


def save_as_alignment(fastadict, path):
    fd = open(path, "wt")
    for seq in fastadict.keys():
        fd.write(seq + "\t" + fastadict.get(seq) + "\t" + seq + "\n")
    fd.close()