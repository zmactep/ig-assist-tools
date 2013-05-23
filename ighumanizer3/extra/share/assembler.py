__author__ = 'mactep'

from Bio.Seq import Seq

nucleotides = "ACGT"


class Vertex(object):
    def __init__(self):
        self.edges = []
        self.in_degree = 0
        self.out_degree = 0


class Edge(object):
    def __init__(self, seq, kmer1, kmer2):
        self.seq = seq
        self.kmer1 = kmer1
        self.kmer2 = kmer2

    def __len__(self):
        return len(self.seq)


def get_kmer_hash(seqs, k, no_reverse=False):
    kmer_hash = {}
    for seq in seqs:
        n = len(seq)
        for i in range(n - k + 1):
            kmer = seq[i:i+k]
            kmer_hash[kmer] = kmer_hash.get(kmer, 0) + 1
            if no_reverse:
                kmer = str(Seq(kmer).reverse_complement())
                kmer_hash[kmer] = kmer_hash.get(kmer, 0) + 1
    return kmer_hash


def get_kmer_next(kmer, kmer_hash):
    return [(kmer[1:] + n) for n in nucleotides if (kmer[1:] + n) in kmer_hash]


def get_kmer_prev(kmer, kmer_hash):
    return [(n + kmer[:-1]) for n in nucleotides if (n + kmer[:-1]) in kmer_hash]


def go_left(kmer, kmer_hash, visited_hash):
    left_edge = []
    prev_list = get_kmer_prev(kmer, kmer_hash)
    while len(prev_list) == 1:
        if prev_list[0] == kmer:
            left_edge = "".join(reversed(left_edge))
            return False, left_edge
        left_edge.append(prev_list[0][0])
        visited_hash[prev_list[0]] = True
        prev_list = get_kmer_prev(prev_list[0], kmer_hash)
    left_edge = "".join(reversed(left_edge))
    return left_edge[:len(kmer) - 1], left_edge[len(kmer) - 1:]


def go_right(kmer, kmer_hash, visited_hash):
    right_edge = []
    next_list = get_kmer_next(kmer, kmer_hash)
    while len(next_list) == 1:
        if next_list[0] == kmer:
            right_edge = "".join(right_edge)
            return False, right_edge
        right_edge.append(next_list[0][-1])
        visited_hash[next_list[0]] = True
        next_list = get_kmer_next(next_list[0], kmer_hash)
    right_edge = "".join(right_edge)
    return right_edge[-(len(kmer) - 1):], right_edge[:-(len(kmer) - 1)]


def buildGraph(k, no_reverse, *seqs):
    kmer_hash = get_kmer_hash(seqs, k + 1, no_reverse)
    print(len(kmer_hash))
    visited_hash = {i: False for i in kmer_hash}
    graph = {}
    for kmer in kmer_hash:
        if visited_hash[kmer]:
            continue
        next_list = get_kmer_next(kmer, kmer_hash)
        prev_list = get_kmer_prev(kmer, kmer_hash)
        if len(next_list) * len(prev_list) == 1:
            # TODO: make this code work
            pass
            # edge = Edge(left_edge + kmer + right_edge, left_vertex, right_vertex)
            # graph[left_vertex].edges.append(edge)
            # graph[left_vertex].out_degree += 1
            # graph[right_vertex].in_degree += 1
        else:
            if len(prev_list) != 1 and kmer[:-1] not in graph:
                graph[kmer[:-1]] = Vertex()
            if len(next_list) != 1 and kmer[1:] not in graph:
                graph[kmer[1:]] = Vertex()
        visited_hash[kmer] = True
    return graph


def test():
    from ighumanizer3.extra.share import fasta_tools
    h3b = fasta_tools.read_abi("/home/mactep/BIO/Ig/abi/E3dMP1_10_H3b_2013-04-17_EGFR_E3dMP1-H3b_002_1.ab1")
    seqr = fasta_tools.read_abi("/home/mactep/BIO/Ig/abi/E3dMP1_10_SeqR_2013-04-17_EGFR_E3dMP1-SeqR_002_A10.ab1")
    s_h3b = h3b.get(list(h3b.keys())[0])
    s_seqr = seqr.get(list(seqr.keys())[0])

    g = buildGraph(11, s_h3b, s_seqr)
    print(len(g))
