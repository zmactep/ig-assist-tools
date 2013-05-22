#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = 'mactep'

from Bio import SeqIO

########################################################


class FastaDict(object):
    def __init__(self):
        self.seq_dict = {}

    def get(self, seq_name):
        return self.seq_dict[seq_name]

    def set(self, seq_name, seq_value):
        self.seq_dict[seq_name] = seq_value

    def contains(self, seq_name):
        return seq_name in self.seq_dict

    def keys(self):
        return self.seq_dict.keys()


class FastqDict(FastaDict):
    def __init__(self):
        super(FastaDict, self).__init__()
        self.qual_dict = {}

    def set(self, seq_name, seq_value):
        self.seq_dict[seq_name] = seq_value
        self.qual_dict[seq_name] = []

    def setqual(self, seq_name, qual_value):
        if self.contains(seq_name):
            self.qual_dict[seq_name] = qual_value

    def setq(self, seq_name, seq_value, qual_value):
        self.set(seq_name, seq_value)
        self.setqual(seq_name, qual_value)

    def getqual(self, seq_name):
        return self.qual_dict[seq_name]


########################################################


def read_fasta(filename, filter_gaps=True):
    result = FastaDict()
    data = SeqIO.parse(filename, "fasta")
    for seq in data:
        s = str(seq.seq)
        if filter_gaps:
            s = s.replace('-', '')
        result.set(seq.id, s)
    return result


def read_qual(filename, qualtype, filter_gaps=True):
    result = FastqDict()
    data = SeqIO.parse(filename, qualtype)
    for seq in data:
        s = str(seq.seq)
        if filter_gaps:
            s = s.replace('-', '')
        result.setq(seq.id, s, seq.letter_annotations["phred_quality"])
    return result


def read_fastq(filename, filter_gaps=True):
    return read_qual(filename, "fastq", filter_gaps)


def read_abi(filename, filter_gaps=True):
    return read_qual(filename, "abi", filter_gaps)

########################################################


def write_fasta(filename, data, line_length=60):
    if isinstance(filename, str):
        fd = open(filename, "w")
    else:
        fd = filename
    for i in data.keys():
        fd.write(">" + i + "\n")
        s = 0
        while s < len(data.get(i)):
            fd.write(data.get(i)[s:s + line_length] + "\n")
            s += line_length
        fd.write("\n")
    if isinstance(filename, str):
        fd.close()


def write_fastq(filename, data):
    if isinstance(filename, str):
        fd = open(filename, "w")
    else:
        fd = filename
    for i in data.keys():
        fd.write("@" + i + "\n")
        fd.write(data.get(i) + "\n")
        fd.write("+\n")
        for q in data.getqual(i):
            fd.write(chr(q))
        fd.write("\n")
    if isinstance(filename, str):
        fd.close()


########################################################


def split_fast(filename, ftype, file_length):
    fd = open(filename, "r")
    data = SeqIO.parse(fd, ftype)
    basename = filename[:filename.rfind('.')]
    seq_list = []
    counter = 0
    for i in data:
        seq_list.append(i)
        if len(seq_list) == file_length:
            counter += 1
            fdo = open(basename + str(counter) + "." + ftype, "w")
            SeqIO.write(seq_list, fdo, ftype)
            fdo.close()
    if seq_list:
        counter += 1
        fdo = open(basename + str(counter) + "." + ftype, "w")
        SeqIO.write(seq_list, fdo, ftype)
        fdo.close()


def split_fasta(filename, file_length):
    split_fast(filename, "fasta", file_length)


def split_fastq(filename, file_length):
    split_fast(filename, "fastq", file_length)