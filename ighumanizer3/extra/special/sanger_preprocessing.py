#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = 'mactep'

from ighumanizer3.extra.share import fasta_tools


def read_and_filter(filename):
    s = fasta_tools.read_abi(filename)
    seq = s.get(s.keys()[0])
    return seq