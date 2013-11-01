# -*- coding: utf-8 -*-

__author__ = 'mactep'

import os
import logging
from Bio.Seq import Seq
from ighumanizer3.extra.share import fasta_tools, algorithm


VL_LEADER_DEFAULT = "ATGAAATACCTGCTGCCGACCGCTGCTGCTGGTCTGCTGCTCCTCGCTGCCCAGCCGGCGATGGCTAGC"
VH_LEADER_DEFAULT = "ATGAAATACCTATTGCCTACGGCAGCCGCTGGATTGTTATTACTCGCGGCCCAGCCGGCCATGGCC"

SANGER_ERROR = "NNNNN"


def get_basename(filename_s, filename_t):
    basename = []
    bn_s = os.path.basename(filename_s)
    bn_t = os.path.basename(filename_t)
    i = 0
    while bn_s[i] == bn_t[i]:
        basename.append(bn_s[i])
        i += 1
    return "".join(basename)


def get_igv(sequence, vl_leader, vh_leader):
    vl_start = sequence.find(vl_leader) + len(vl_leader)
    vl_end = sequence.find(vh_leader)
    vh_start = vl_end + len(vh_leader)
    return sequence[vl_start:vl_end], sequence[vh_start:]


def get_igva(vl, vh):
    vla = str(Seq(vl).translate())
    vha = str(Seq(vh).translate())
    return vla[:vla.find('*')], vha[:vha.find('*')]


def construct_strings(sanger_first, sanger_second, reverse_complement=False):
    s = sanger_first.get(list(sanger_first.keys())[0])
    t = sanger_second.get(list(sanger_second.keys())[0])
    if reverse_complement:
        t = str(Seq(t).reverse_complement())
    return s, t


def construct_quality(sanger_first, sanger_second, reverse_complement=False):
    s = sanger_first.getqual(list(sanger_first.keys())[0])
    t = sanger_second.getqual(list(sanger_second.keys())[0])
    if reverse_complement:
        t = reversed(t)
    return s, t


def merge_naively(s, t, s_qual=None, t_qual=None):
    res = algorithm.longest_common_substring(s, t)
    if res:
        (s_start, s_end), (t_start, t_end) = res
        return s[:s_end] + t[t_end:]
    return False


def read_and_merge(filename_s, filename_t, merge_function):
    logging.debug("read_and_merge")
    logging.debug("Reading %s" % filename_s)
    sanger_s = fasta_tools.read_abi(filename_s)
    logging.debug("Sequence %s readed" % filename_s)
    logging.debug("Reading %s" % filename_t)
    sanger_t = fasta_tools.read_abi(filename_t)
    logging.debug("Sequence %s readed" % filename_t)
    s, t = construct_strings(sanger_s, sanger_t, True)
    s_qual, t_qual = construct_quality(sanger_s, sanger_t, True)
    logging.debug("S: %s" % s)
    logging.debug("T: %s" % t)
    if s == SANGER_ERROR and t != SANGER_ERROR:
        return t
    elif s != SANGER_ERROR and t == SANGER_ERROR:
        return s
    elif s == SANGER_ERROR and t == SANGER_ERROR:
        return False
    else:
        return merge_function(s, t, s_qual, t_qual)


def read_and_process(filename_s, filename_t, forward_mark, backward_mark, vl_leader, vh_leader,
                     merge_function=merge_naively):
    logging.debug("read_and_process")
    sequence = None
    if forward_mark in filename_s and backward_mark in filename_t:
        sequence = read_and_merge(filename_s, filename_t, merge_function)
    elif forward_mark in filename_t and backward_mark in filename_s:
        sequence = read_and_merge(filename_t, filename_s, merge_function)
    if sequence and vl_leader in sequence and vh_leader in sequence:
        vl, vh = get_igv(sequence, vl_leader, vh_leader)
        vla, vha = get_igva(vl, vh)
        basename = get_basename(filename_s, filename_t)
        fd_n = fasta_tools.FastaDict()
        fd_n.set(basename + "nucleo-VL", vl)
        fd_n.set(basename + "nucleo-VH", vh)
        fd_a = fasta_tools.FastaDict()
        fd_a.set(basename + "amino-VL", vla)
        fd_a.set(basename + "amino-VH", vha)
        print("Ok")
        return fd_n, fd_a
    else:
        print("Fail")
        return False


def read_and_write(filename_s, filename_t, forward_mark, backward_mark, vl_leader, vh_leader, out_dir,
                   merge_function=merge_naively):
    logging.debug("read_and_write")
    fd = read_and_process(filename_s, filename_t, forward_mark, backward_mark, vl_leader, vh_leader, merge_function)
    if fd:
        basename = get_basename(filename_s, filename_t)
        path_n = os.path.join(out_dir, basename + "-nucleo.fa")
        path_a = os.path.join(out_dir, basename + "-amino.fa")
        fasta_tools.write_fasta(path_n, fd[0])
        fasta_tools.write_fasta(path_a, fd[1])


def run_on_directory(directory, out_directory, fm="SeqR", bm="H3b"):
    logging.debug("run_on_directory")
    cache = []
    filenames = list(filter(lambda f: f.endswith(".ab1"), os.listdir(directory)))
    filenames.sort()
    for filename in filenames:
        s = os.path.join(directory, filename)
        cache.append(s)
        logging.debug("Added to cache: %s" % s)
        if len(cache) == 2:
            logging.debug("Cache: {0} {1}".format(cache[0], cache[1]))
            if cache[0][:cache[0].find(fm)] != cache[1][:cache[1].find(bm)] and \
               cache[0][:cache[0].find(bm)] != cache[1][:cache[1].find(fm)]:
                logging.debug("Cache names error. Shifting cache.")
                cache = [cache[1]]
                continue
            logging.debug("Processing: %s %s" % (cache[0], cache[1]))
            try:
                read_and_write(cache[0], cache[1], fm, bm, VL_LEADER_DEFAULT, VH_LEADER_DEFAULT, out_directory)
            except:
                logging.debug("Bad files, skipping.")
                print("Fail")
                cache = []
                continue
            cache = []