# -*- coding: utf-8 -*-

__author__ = 'mactep'

import os
import logging
from ighumanizer3.extra.share.fasta_tools import FastaDict, read_fasta
from ighumanizer3.extra.share.alignment import *

CHAIN_MARKER_VL = "VL"
CHAIN_MARKER_VH = "VH"


def load_chains(directory, suffix):
    filenames = list(filter(lambda f: f.endswith(suffix), os.listdir(directory)))
    heavy_chains = FastaDict()
    light_chains = FastaDict()
    logging.debug("Creating VL and VH arrays")
    for filename in filenames:
        name = os.path.join(directory, filename)
        fd = read_fasta(name)
        for seq in fd.keys():
            if seq.endswith(CHAIN_MARKER_VL):
                light_chains.set(seq, fd.get(seq))
            elif seq.endswith(CHAIN_MARKER_VH):
                heavy_chains.set(seq, fd.get(seq))
    return light_chains, heavy_chains


def load_directory(directory):
    light_chains_a, heavy_chains_a = load_chains(directory, "-amino.fa")
    light_chains_n, heavy_chains_n = load_chains(directory, "-nucleo.fa")

    logging.debug("Aligning VL")
    light_chains = multiple_alignment(light_chains_a, SeqTypeData().TYPE_UNI_FAST)
    logging.debug("Aligning VH")
    heavy_chains = multiple_alignment(heavy_chains_a, SeqTypeData().TYPE_UNI_FAST)
    logging.debug("Alignment ok!")

    return light_chains_a, heavy_chains_a, light_chains_n, heavy_chains_n


def save_as_alignment(fastadict, path):
    fd = open(path, "wt")
    for seq in fastadict.keys():
        fd.write(seq + "\t" + fastadict.get(seq) + "\t" + seq + "\n")
    fd.close()