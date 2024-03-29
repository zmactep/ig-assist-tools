# -*- coding: utf-8 -*-

__author__ = 'mactep'

import os
import sys
import logging
from ighumanizer3.extra.share.fasta_tools import FastaDict, read_fasta, write_fasta
from ighumanizer3.extra.share.alignment import *
from ighumanizer3.extra.share.algorithm import get_common_name

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


def load_directory(directory, file_usage=False):
    light_chains_a, heavy_chains_a = load_chains(directory, "-amino.fa")
    light_chains_n, heavy_chains_n = load_chains(directory, "-nucleo.fa")

    common = get_common_name(directory, "-amino.fa")

    write_fasta(os.path.join(directory, common + "light-chains-amino.fa"), light_chains_a)
    write_fasta(os.path.join(directory, common + "heavy-chains-amino.fa"), heavy_chains_a)
    write_fasta(os.path.join(directory, common + "light-chains-nucleo.fa"), light_chains_n)
    write_fasta(os.path.join(directory, common + "heavy-chains-nucleo.fa"), heavy_chains_n)

    logging.debug("Aligning VL")
    if sys.platform != "win32" and not file_usage:
        light_chains_a = multiple_alignment(light_chains_a, SeqTypeData().TYPE_UNI_FAST)
        light_chains_n = multiple_alignment(light_chains_n, SeqTypeData().TYPE_UNI_FAST)
    else:
        light_chains_a = multiple_alignment_use_files(os.path.join(directory, common + "light-chains-amino.fa"), SeqTypeData().TYPE_UNI_FAST)
        light_chains_n = multiple_alignment_use_files(os.path.join(directory, common + "light-chains-nucleo.fa"), SeqTypeData().TYPE_UNI_FAST)
    logging.debug("Aligning VH")
    if sys.platform != "win32" and not file_usage:
        heavy_chains_a = multiple_alignment(heavy_chains_a, SeqTypeData().TYPE_UNI_FAST)
        heavy_chains_n = multiple_alignment(heavy_chains_n, SeqTypeData().TYPE_UNI_FAST)
    else:
        heavy_chains_a = multiple_alignment_use_files(os.path.join(directory, common + "heavy-chains-amino.fa"), SeqTypeData().TYPE_UNI_FAST)
        heavy_chains_n = multiple_alignment_use_files(os.path.join(directory, common + "heavy-chains-nucleo.fa"), SeqTypeData().TYPE_UNI_FAST)
    logging.debug("Alignment ok!")

    return light_chains_a, heavy_chains_a, light_chains_n, heavy_chains_n, common


def save_as_alignment(fastadict, path):
    fd = open(path, "wt")
    for seq in fastadict.keys():
        fd.write(seq + "\t" + fastadict.get(seq) + "\t" + seq + "\n")
    fd.close()