#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = 'mactep'

# Matrices
from Bio.SubsMat import MatrixInfo


class Matrices(object):
    def __init__(self):
        self.blosum62 = MatrixInfo.blosum62
        self.pam250 = MatrixInfo.pam250

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

import os
import sys
import copy
import subprocess

from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from ighumanizer3.extra.share import fasta_tools
if sys.version_info[0] == 3:
    from io import StringIO
    import iterpipes3 as iterpipes
else:
    from StringIO import StringIO
    import iterpipes


class SeqTypeData(object):
    def __init__(self):
        self.TYPE_DEFAULT = 0
        self.TYPE_UNI_FAST = 1
        self.TYPE_AMINO_FAST = 2
        self.TYPE_NUCLEO_FAST = 3
        self.type2cmd = {self.TYPE_DEFAULT: MuscleCommandline(clwstrict=True),
                         self.TYPE_UNI_FAST: MuscleCommandline(clwstrict=True, maxiters=2),
                         self.TYPE_AMINO_FAST: MuscleCommandline(clwstrict=True, maxiters=1,
                                                                 diags=True, sv=True, distance1="kbit20_3"),
                         self.TYPE_NUCLEO_FAST: MuscleCommandline(clwstrict=True, maxiters=1, diags=True)}
        self.type2cmd_f = {self.TYPE_DEFAULT: MuscleCommandline(),
                           self.TYPE_UNI_FAST: MuscleCommandline(maxiters=2),
                           self.TYPE_AMINO_FAST: MuscleCommandline(maxiters=1,
                                                                   diags=True, sv=True, distance1="kbit20_3"),
                           self.TYPE_NUCLEO_FAST: MuscleCommandline(maxiters=1, diags=True)}


def multiple_alignment_use_files(file_input, alignment_type=SeqTypeData().TYPE_DEFAULT):
    muscle_cmd = SeqTypeData().type2cmd_f[alignment_type]
    result = iterpipes.run(iterpipes.linecmd(str(muscle_cmd) + " -in \"" + file_input + "\""))
    result_list = []
    for i in result:
        result_list.append(i.strip())
    # iores = StringIO()
    # for i in result_list:
    #     print(i[0])
    #     iores.write(i)
    name = file_input[:file_input.rfind('.')] + "-aligned.fa"
    with open(name, "wt") as tmpfd:
        tmpfd.write("\n".join(result_list))
    return fasta_tools.read_fasta(name, False)


def multiple_alignment(fasta_dict, alignment_type=SeqTypeData().TYPE_DEFAULT):
    in_handle = StringIO()
    fasta_tools.write_fasta_handle(in_handle, fasta_dict)

    muscle_cmd = SeqTypeData().type2cmd[alignment_type]
    child = subprocess.Popen(str(muscle_cmd), stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             shell=(sys.platform != "win32"))
    if not child:
        print("Process was not created!")
        return

    if sys.version_info[0] == 3:
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