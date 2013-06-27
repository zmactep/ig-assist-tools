# -*- coding: utf-8 -*-

__author__ = 'mactep'

import os
import data.train
import numpy as np
from itertools import chain
from sklearn import svm
from ighumanizer3.extra.svm.encodimg_tools import encode_area


class RegionsData(object):
    def __init__(self):
        self.num2region = [i + str(j) for j in range(1, 4) for i in ["FR", "CDR"]] + ["FR4", "C"]
        self.train_set = os.path.join(os.path.abspath(data.train.__path__[0]), "regions")


class TrainRegionsDataset(object):
    def __init__(self, radius=8):
        self.sequences = []
        self.regions = []
        self.radius = radius

    def add(self, sequences, regions):
        self.sequences.append(sequences)
        self.regions.append(regions)

    def merge(self, train_dataset, save_regions=True):
        for sequences in train_dataset.sequences:
            self.sequences.append(sequences)
        for regions in train_dataset.regions:
            self.regions.append(regions)
        if not save_regions:
            self.radius = train_dataset.radius

    def next(self):
        for pos, sequences in enumerate(self.sequences):
            regions = self.regions[pos]
            for seq in sequences:
                s = seq.replace(' ', '').replace('-', '')
                for i, c in enumerate(seq):
                    if c != ' ' and c != '':
                        delta = seq[:i].count(' ') + seq[:i].count('-')
                        yield encode_area(s, i - delta, self.radius), int(regions[i])-1


class RegionsClassifier(object):
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

    def train(self, path):
        d = filter(lambda s: s.endswith(".train"), os.listdir(path))
        dataset = TrainRegionsDataset(self.radius)
        for fn in d:
            filename = os.path.join(path, fn)
            dataset.merge(construct_train_regions_set(filename, self.radius))
        _input = []
        _output = []
        for _in, _out in dataset.next():
            _input.append(_in)
            _output.append(_out)
        self.classifier.fit(np.asarray(list(chain(*_input)), float), np.asarray(_output, float))

    def predict(self, protein):
        result = []
        for i, _ in enumerate(protein):
            result.append(list(self.classifier.predict(encode_area(protein, i, self.radius))))
        return "".join(str(int(i)) for i in chain(*result))


def construct_train_regions_set(filename, radius=8):
    sequences = []
    fd = open(filename, "rt")
    regions = fd.readline().strip()
    for line in fd:
        sequences.append(line.strip())
    fd.close()
    t = TrainRegionsDataset(radius)
    t.add(sequences, regions)
    return t


def euristic_fix(t):
    seq = []
    for i, c in enumerate(t):
        ic = int(c)
        if i:
            ip = int(t[i-1])
            if ic != ip + 1 and ic != ip:
                ic = ip
        seq.append(str(ic))
    return "".join(seq)


def test(n=None):
    seq = ["DVVMTQSPSSVTASVGETVTISCKSSQSVAYKSNQKNYLAWYQQRPGQSPRLLIYWASTRTPGIPDRFSGNGSTTDFTMTISSFQPEDAA",
           "QSVLTQPPSVSGTLGNTVTIACAGTSSDVGSGNYVSWYQQLPGTAPKTIIYQDNKRLPGIPDRFSGSKSGNTAFLTISGLQSLDDADYYC",
           "QLVLTQSPSASASLGTAVKLTCTLNSQYSTYYIHWFQQKLGQTPSFLMKVTSDGRVVKGDGVPGRFSGSSSGADRYLTVSNIQSEDEADY",
           "VQSGGGLVQAGGSLRLSCAASGGTFATSPMGWLRQAPGKEREFVAAISPSGGDRIYDDSVKGRFTISRDNAGYFIYLQMNSLKPEDTARY",
           "DVVMTQSPSSVTASAGETASINCKSSQSVLYSSNQKNYLAWYQQRPGQSPRLLIYWASTRESGVPDRFSGSGSTTDFTLTISSFQPEDAA",
           "DIQMTQSPASLSHSLGDRVDITCQASQSINNKIAWYQQKPGHPPKVLIYAASKLPTGVPSRFSGSGSGTTFTLTINELEAQDVATYYCLQ",
           "QAVLTQPPSVSGSPGQTVTISCTGTSDDVGSGNYVSWYQQVPGMAPKLLIYNAGTRRAGITGRFSASKSGNTASLTISGLQSEDEADYYC"]
    tru = ["000000000000000000000001111111111111111122222222222222233333334444444444444444444444444444",
           "000000000000000000000011111111111111222222222222222333333344444444444444444444444444444555",
           "000000000000000000000011111111111122222222222222233333333333444444444444444444444444444445",
           "000000000000000000000000001111122222222222222333333333333333334444444444444444444444444445",
           "000000000000000000000001111111111111111122222222222222233333334444444444444444444444444444",
           "000000000000000000000001111111111122222222222222233333334444444444444444444444444444455555",
           "000000000000000000000011111111111111222222222222222333333344444444444444444444444444444555"]
    result = []
    align = []

    def run(i):
        print("Testing window radius %i (width: %i)" % (i + 1, 2 * (i + 1) + 1))
        r = RegionsClassifier(i + 1)
        r.train(RegionsData().train_set)
        result.append([])
        align.append([])
        for j, s in enumerate(seq):
            k = tru[j]
            t = r.predict(s)
            c = sum(int(k[i] == t[i]) for i, _ in enumerate(t)) / len(t)
            result[-1].append(c)
            align[-1].append((k, t))

    if not n:
        for i in range(20):
            run(i)
    else:
        run(n)

    return result, align