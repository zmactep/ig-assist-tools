__author__ = 'mactep'

import os
import random
import data.train
import numpy as np
from Bio import SeqIO
from collections import Counter
from itertools import chain
from sklearn import svm
from ighumanizer3.extra.svm.encodimg_tools import encode_area


class TypeData(object):
    def __init__(self):
        self.num2type = ["VL", "VK", "VH", "HH"]
        self.train_set = os.path.join(os.path.abspath(data.train.__path__[0]), "type")


class TrainTypeDataset(object):
    def __init__(self, radius):
        self.types = [[] for _ in TypeData().num2type]
        self.radius = radius

    def add(self, sequences, tp):
        for seq in sequences:
            self.types[tp].append(seq)

    def merge(self, train_dataset):
        for i, tp in enumerate(train_dataset.types):
            for seq in tp:
                self.types[i].append(seq)

    def next(self):
        for i, tp in enumerate(self.types):
            for seq in tp:
                s = seq.replace(' ', '').replace('-', '')
                for pos, c in enumerate(s):
                    yield encode_area(s, pos, self.radius), i


class TypeClassifier(object):
    def __init__(self, radius):
        self.classifier = svm.SVC(kernel='rbf')
        self.radius = radius

    def train_dataset(self, dataset):
        _input = []
        _output = []
        for _in, _out in dataset.next():
            _input.append(_in)
            _output.append(_out)
        self.classifier.fit(np.asarray(list(chain(*_input)), float), np.asarray(_output, float))

    def train(self, path, fasta=False):
        d = filter(lambda s: s.endswith(".train" + ".fa" * int(fasta)), os.listdir(path))
        dataset = TrainTypeDataset(self.radius)
        for fn in d:
            tp = TypeData().num2type.index(fn[:2].upper())
            filename = os.path.join(path, fn)
            if not fasta:
                dataset.merge(construct_train_type_set(filename, tp, self.radius))
            else:
                dataset.merge(construct_train_type_set_fasta(filename, tp, self.radius))
        _input = []
        _output = []
        for _in, _out in dataset.next():
            _input.append(_in)
            _output.append(_out)
        self.classifier.fit(np.asarray(list(chain(*_input)), float), np.asarray(_output, float))

    def predict(self, protein):
        result = []
        for i, c in enumerate(protein):
            result.append(list(self.classifier.predict(encode_area(protein, i, self.radius))))
        result = [int(i) for i in chain(*result)]
        return Counter(result).most_common(1)[0][0]


def construct_train_type_set(filename, tp, radius):
    sequences = []
    fd = open(filename, "rt")
    for line in fd:
        sequences.append(line.strip())
    fd.close()
    t = TrainTypeDataset(radius)
    t.add(sequences, tp)
    return t


def construct_train_type_set_fasta(filename, tp, radius):
    sequences = []
    for s in SeqIO.parse(filename, "fasta"):
        sequences.append(str(s.seq))
    t = TrainTypeDataset(radius)
    t.add(sequences, tp)
    return t


def test(n=8):
    result = []
    path = TypeData().train_set
    print("Testing window radius %i (width: %i)" % (n, 2 * n + 1))
    t = TypeClassifier(n)
    print("Training SVM")
    t.train(path, True)
    for tp in ["vh", "vl", "vk"]:
        print("Validating %s results" % tp)
        l = []
        for j in SeqIO.parse(os.path.join(os.path.join(path, "big"), tp + ".train.fa"), "fasta"):
            if random.random() < 0.1:
                l.append(str(j.seq))
        lr = {}
        for j in l:
            pr = TypeData().num2type[t.predict(j)]
            lr[pr] = lr.get(pr, 0) + 1
        result.append(lr)
    return result