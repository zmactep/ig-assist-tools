#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = 'mactep'

# Matrices
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
if sys.version_info.major == 3:
    from io import StringIO
else:
    from StringIO import StringIO

def multiple_alignment(fasta_dict):
    in_handle = StringIO()
    fasta_tools.write_fasta_handle(in_handle, fasta_dict)

    muscle_cmd = MuscleCommandline(clwstrict=True)
    child = subprocess.Popen(str(muscle_cmd), stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             shell=(sys.platform != "win32"))
    if not child:
        print("Process was not created!")
        return

    if sys.version_info.major == 3:
        child.stdin.write(bytes(in_handle.getvalue(), 'utf-8'))
        child.stdin.close()
        align = AlignIO.read(StringIO("".join(line.decode() for line in child.stdout)), "clustal")
    else:
        child.stdin.write(in_handle.getvalue())
        child.stdin.close()
        align = AlignIO.read(child.stdout, "clustal")
    fd = copy.deepcopy(fasta_dict)
    for a in align:
        fd.set(a.id, str(a.seq))

    return fd