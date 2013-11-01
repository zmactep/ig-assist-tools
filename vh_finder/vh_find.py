__author__ = 'mactep'

import sys
import os

from Bio.Seq import Seq
from Bio import SeqIO

from ighumanizer3.extra.share.alignment import *
from ighumanizer3.extra.share.algorithm import get_common_name
from sanger_processing.sanger_processing import VH_LEADER_DEFAULT, VL_LEADER_DEFAULT
from sanger_processing.cluster_analysis import save_as_alignment


def try_hl(rec, seq, leader):
    pos = seq.find(leader)
    if pos != -1:
        return SeqIO.SeqRecord(Seq(seq[pos+len(leader):]), rec.id, rec.name, rec.description)


def try_v(rec, leader):
    seq = str(rec.seq)
    r = try_hl(rec, seq, leader)
    if r:
        return r
    seq = str(Seq(seq).reverse_complement())
    r = try_hl(rec, seq, leader)
    if r:
        return r


def get_trans(rec):
    seq = str(rec.seq.translate())
    seq = seq[:seq.find('*')]
    return SeqIO.SeqRecord(Seq(seq), rec.id, rec.name, rec.description)


def post_process(directory):
    common = get_common_name(directory, "-amino.fa")

    for d_type in ["heavy", "light"]:
        for a_type in ["amino", "nucleo"]:
            s = "%s%s-chains-%s.fa" % (common, d_type, a_type)
            if not len(list(SeqIO.parse(os.path.join(directory, s), "fasta"))):
                continue
            v = multiple_alignment_use_files(os.path.join(directory, s), SeqTypeData().TYPE_UNI_FAST)
            save_as_alignment(v, os.path.join(directory, s.replace(".fa", "-aligned.txt")))


def process(directory, out_dir, pminlen):
    vh = []
    vl = []
    trash = []
    filenames = list(filter(lambda f: f.endswith(".ab1"), os.listdir(directory)))
    for filename in filenames:
        print(os.path.join(directory, filename))
        try:
            s = list(SeqIO.parse(os.path.join(directory, filename), "abi"))[0]
        except Exception:
            continue
        h = try_v(s, VH_LEADER_DEFAULT)
        if h:
            if len(h) >= pminlen * 3:
                vh.append(h)
                continue
        l = try_v(s, VL_LEADER_DEFAULT)
        if l:
            if len(l) >= pminlen * 3:
                vl.append(l)
                continue
        trash.append(s)

    SeqIO.write(vh, os.path.join(out_dir, "heavy-chains-nucleo.fa"), "fasta")
    SeqIO.write(vl, os.path.join(out_dir, "light-chains-nucleo.fa"), "fasta")
    SeqIO.write(map(get_trans, vh), os.path.join(out_dir, "heavy-chains-amino.fa"), "fasta")
    SeqIO.write(map(get_trans, vl), os.path.join(out_dir, "light-chains-amino.fa"), "fasta")
    SeqIO.write(trash, os.path.join(out_dir, "trash.fa"), "fasta")

    post_process(out_dir)