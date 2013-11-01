#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = 'mactep'

from collections import Counter

from Bio import SeqIO
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq

########################################################


class FastaDict(object):
    def __init__(self):
        self.seq_dict = {}

    def get(self, seq_name):
        return self.seq_dict[seq_name]

    def set(self, seq_name, seq_value):
        self.seq_dict[seq_name] = str(seq_value)

    def contains(self, seq_name):
        return seq_name in self.seq_dict

    def keys(self):
        return self.seq_dict.keys()

    def __contains__(self, item):
        return item in self.seq_dict

    def __getitem__(self, item):
        return self.seq_dict[item]

    def __setitem__(self, key, value):
        self.seq_dict[key] = value


class FastqDict(FastaDict):
    def __init__(self):
        super(FastaDict, self).__init__()
        self.seq_dict = {}
        self.qual_dict = {}

    def set(self, seq_name, seq_value):
        self.seq_dict[seq_name] = str(seq_value)
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


def get_consensus_letter(fasta, xletter):
    result = []
    keys = fasta.keys()
    ll = len(fasta.get(keys[0]))
    for key in keys:
        if len(fasta.get(key)) != ll:
            return ""
    for i in range(ll):
        letter, count = Counter(fasta.get(key)[i] for key in keys).most_common(1)[0]
        if count >= len(keys) * 0.5:
            result.append(letter)
        else:
            result.append(xletter)
    return "".join(result)


def get_consensus_amino(fasta):
    return get_consensus_letter(fasta, 'X')


def get_consensus_nucleo(fasta):
    return get_consensus_letter(fasta, 'N')


########################################################


def fasta_from_keylist(fasta, keylist):
    return {key: fasta.get(key) for key in keylist if key in fasta}


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


def write_fasta(filename, data):
    fd = open(filename, "w")
    seq_list = []
    for i in data.keys():
        seq_list.append(SeqRecord(Seq(data.get(i)), id=i, description=""))
    SeqIO.write(seq_list, fd, "fasta")
    fd.close()


def write_fastq(filename, data):
    fd = open(filename, "w")
    seq_list = []
    for i in data.keys():
        seq_list.append(SeqRecord(Seq(data.get(i)), id=i, description="",
                        letter_annotations={'solexa_quality': data.getqual(i)}))
    SeqIO.write(seq_list, fd, "fastq")
    fd.close()


########################################################


def write_fasta_handle(handle, data):
    seq_list = []
    for i in data.keys():
        seq_list.append(SeqRecord(Seq(data.get(i)), id=i, description=""))
    SeqIO.write(seq_list, handle, "fasta")
    

def write_fastq_handle(handle, data):
    seq_list = []
    for i in data.keys():
        seq_list.append(SeqRecord(Seq(data.get(i)), id=i, description="",
                        letter_annotations={'solexa_quality': data.getqual(i)}))
    SeqIO.write(seq_list, handle, "fastq")

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