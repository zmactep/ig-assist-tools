#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = 'mactep'

# Matrices
from sys import version_info

if version_info.major == 3:
    from io import StringIO

from Bio.SubsMat import MatrixInfo

blosum62 = MatrixInfo.blosum62
pam250 = MatrixInfo.pam250

#####################################################

# Pairwise

from Bio import pairwise2


def global_alignment(s, t, matrix, gap_open, gap_ext):
    return [(s_gaps, t_gaps, score)
            for s_gaps, t_gaps, score, start, end
            in pairwise2.align.globalds(s, t, matrix, gap_open, gap_ext)]


def local_alignment(s, t, matrix, gap_open, gap_ext):
    return pairwise2.align.localds(s, t, matrix, gap_open, gap_ext)


def edit_distance(s, t):
    return [(s_gaps, t_gaps, score) for s_gaps, t_gaps, score, start, end in pairwise2.align.globalxx(s, t)]

#####################################################

# Multiple

import sys
import copy
import subprocess

from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from ighumanizer3.extra.share import fasta_tools


def multiple_alignment(fasta_dict):
    # fasta_tools.write_fasta_handle(sys.stdout, fasta_dict)
    # print("******************************")
    muscle_cmd = MuscleCommandline(clwstrict=True)
    child = subprocess.Popen(str(muscle_cmd), stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             shell=(sys.platform != "win32"))
    if not child:
        print("Process was not created!")
        return

    # print("Writing data to MUSCLE")
    fasta_tools.write_fasta_handle(child.stdin, fasta_dict)
    child.stdin.close()
    # print("Data was written")

    if version_info.major == 3:
        align = AlignIO.read(StringIO("".join(line.decode() for line in child.stdout)), "clustal")
    else:
        align = AlignIO.read(child.stdout, "clustal")
    fd = copy.deepcopy(fasta_dict)
    for a in align:
        fd.set(a.id, str(a.seq))

    return fd