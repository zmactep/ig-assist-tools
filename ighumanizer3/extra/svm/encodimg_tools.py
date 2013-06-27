# -*- coding: utf-8 -*-

__author__ = 'mactep'

import itertools
from Bio.Data.CodonTable import unambiguous_dna_by_id


def codon2bits(codon):
    r = []
    for num in codon:
        r.append([int(num == i) for i in [1, 2, 4, 8]])
    return list(itertools.chain(*r))


class Encoder(object):
    multiflag = True

    def __init__(self):
        tmp = unambiguous_dna_by_id[1].forward_table
        self.bits = {'A': 1, 'C': 2, 'G': 4, 'T': 8}
        self.table_multi = {}
        for k in tmp:
            l = self.table_multi.get(tmp[k], [])
            l.append(list(map(self.bits.__getitem__, k)))
            self.table_multi[tmp[k]] = l
        self.table_multi['*'] = [[0, 0, 0]]
        tmp = unambiguous_dna_by_id[1].back_table
        self.table_single = {k: [list(map(self.bits.__getitem__, tmp[k]))] for k in tmp}
        self.table_single['*'] = [[0, 0, 0]]

    def __new__(cls, *args, **kwargs):
        if not hasattr(cls, 'instance'):
            cls.instance = super(Encoder, cls).__new__(cls)
        return cls.instance

    def __getitem__(self, item):
        result = []
        table = self.table_multi if self.multiflag else self.table_single
        for alpha in item:
            new_result = []
            for codon in table[alpha]:
                e = codon2bits(codon)
                new_result += [r + e for r in result] if result else [e]
            result = new_result
        return result

    def get(self, item, flag):
        state = self.multiflag
        self.multiflag = flag
        x = self[item]
        self.multiflag = state
        return x

    def get_multi(self, item):
        return self.get(item, True)

    def get_single(self, item):
        return self.get(item, False)


def encode_protein(protein, multi=False):
    a = Encoder()
    return a.get(protein, multi)


def encode_area(protein, pos, radius=8, multi=False):
    s = '*' * (radius - pos) + \
        protein[pos-radius if pos-radius > 0 else 0:1+pos+radius] +\
        '*' * (pos + radius + 1 - len(protein))
    return encode_protein(s, multi)