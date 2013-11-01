__author__ = 'mactep'

import logging

import os
import difflib
import itertools
import shutil
import logging
import pprint
import textwrap
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO, SeqIO, Phylo
from Bio.Seq import Seq
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Align.Applications import ClustalwCommandline
from Bio.Alphabet import IUPAC
from Bio.SeqUtils.CheckSum import seguid
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

from ighumanizer3.extra.share.fasta_tools import write_fasta, read_fasta
from ighumanizer3.extra.share.algorithm import get_common_name

from cluster_processing import style, color_seq


def remove_duplicates(fasta):
    check_sums = set()
    for record in fasta:
        check_sum = seguid(record.seq)
        if check_sum in check_sums:
            logging.info("Ignoring record {0}".format(record.id))
            continue
        check_sums.add(check_sum)
        yield record
        # yield SeqIO.SeqRecord(id="ID=" + record.id,
        #                       description=record.description,
        #                       seq=record.seq, name=record.name)


def distance(record1, record2):
    return 1 - difflib.SequenceMatcher(None, str(record1.seq),
                                       str(record2.seq)).ratio()


def filter_old_gaps(fasta):
    for record in fasta:
        rec = record
        rec.seq = Seq(str(record.seq).replace('-', ''))
        yield rec


def processing(raw_fasta_path, out_dir_path, prct=False):
    if not os.path.exists(out_dir_path):
        logging.info("Making directory {0}".format(out_dir_path))
        os.makedirs(out_dir_path)

    # if prct:
    #     deduplicated_fasta = remove_duplicates(
    #         SeqIO.parse(raw_fasta_path, "fasta"))
    # else:
    #     deduplicated_fasta = SeqIO.parse(raw_fasta_path, "fasta")
    deduplicated_fasta = filter_old_gaps(SeqIO.parse(raw_fasta_path, "fasta"))
    base = os.path.basename(raw_fasta_path)
    fasta_path = os.path.join(out_dir_path, base)

    logging.info("Writing FASTA in {0}".format(fasta_path))
    SeqIO.write(deduplicated_fasta, fasta_path, "fasta")

    # Multiple sequence alignment
    cline = ClustalwCommandline("clustalw2", infile=fasta_path)
    stdout, stderr = cline()
    logging.info(cline)

    if fasta_path.endswith(".fasta"):
        clustalw_result_path = fasta_path.replace(".fasta", ".aln")
    else:
        clustalw_result_path = fasta_path.replace(".fa", ".aln")

    alignment_dict = SeqIO.to_dict(
        AlignIO.read(clustalw_result_path, "clustal"))

    # writing alignment table in .txt
    with open(os.path.join(out_dir_path, "alignment.txt"), "w") as fout:
        fout.write(
            "\n".join(
                str(record.seq) for record in alignment_dict.itervalues()))

    # alignment tree drawing
    if fasta_path.endswith(".fasta"):
        tree_path = fasta_path.replace(".fasta", ".dnd")
    else:
        tree_path = fasta_path.replace(".fa", ".dnd")
    tree = Phylo.read(tree_path, "newick")
    tree.ladderize()

    # with labels
    # Phylo.draw_graphviz(tree, label_func=lambda x: x.name.replace("ID=", ""))
    Phylo.draw_graphviz(tree)
    plt.savefig(os.path.join(out_dir_path,
                               "figure_with_labels.pdf"))  # need pygraphviz, pylab

    # Clustering
    ids = dict(enumerate(alignment_dict.keys()))
    distance_matrix = np.zeros([len(ids)] * 2)
    for i, j in itertools.combinations(xrange(len(ids)), r=2):
        distance_matrix[i][j] = distance_matrix[j][i] = \
            distance(alignment_dict[ids[i]], alignment_dict[ids[j]])

    # Compute and plot dendrogram
    fig = plt.figure()
    axdendro = fig.add_axes([0.09, 0.1, 0.2, 0.8])
    Y = linkage(distance_matrix, method="centroid")
    cutoff = 0.5 * max(Y[:, 2])
    clusters = fcluster(Y, cutoff, "distance")
    Z = dendrogram(Y, orientation="right", color_threshold=cutoff)
    axdendro.set_yticks([])

    # Plot distance matrix
    axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.8])
    index = Z["leaves"]
    distance_matrix = distance_matrix[index, :]
    distance_matrix = distance_matrix[:, index]
    im = axmatrix.matshow(distance_matrix, aspect="auto", origin="lower")
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Plot colorbar
    axcolor = fig.add_axes([0.91, 0.1, 0.02, 0.8])
    plt.colorbar(im, cax=axcolor)

    # Display and save figure
    dendogram_path = os.path.join(out_dir_path, "dendogram.png")
    fig.savefig(dendogram_path)

    fasta_clusters = defaultdict(list)
    for i, cluster in enumerate(clusters):
        fasta_id = ids[i]
        fasta_clusters[cluster].append(alignment_dict[fasta_id])

    # Saving information about clusters
    clusters_dir_path = os.path.join(out_dir_path, "clusters")
    if not os.path.exists(clusters_dir_path):
        os.makedirs(clusters_dir_path)
    freq_dir_path = os.path.join(out_dir_path, "frequencies")
    if not os.path.exists(freq_dir_path):
        os.makedirs(freq_dir_path)
    clusters_meta_path = os.path.join(out_dir_path, "clusters_meta.txt")
    meta_file = open(clusters_meta_path, "w")
    for cluster_id, cluster in fasta_clusters.iteritems():
        cluster_path = os.path.join(clusters_dir_path,
                                    "cluster_{0}.fasta".format(cluster_id))
        SeqIO.write(cluster, cluster_path, "fasta")
        summary_align = AlignInfo.SummaryInfo(MultipleSeqAlignment(cluster))
        consensus = summary_align.dumb_consensus()
        pssm = summary_align.pos_specific_score_matrix(consensus,
                                                       chars_to_ignore=['-'])
        frequencies = dict.fromkeys(IUPAC.protein.letters, 0)
        frequencies.update(
            (key, len(list(group)))
            for key, group in itertools.groupby(sorted(consensus)))
        if "X" in frequencies:
            frequencies.pop("X")

        meta_file.write("""Cluster ID: {0}
Cluster size: {1}
Consensus:
{2}

PSSM:
{3}
Frequencies in consensus:
{4}
""".format(cluster_id, len(cluster), textwrap.fill(str(consensus)), pssm,
           pprint.pformat(frequencies)))

        fig = plt.figure()
        pos = np.arange(len(IUPAC.protein.letters))
        width = .5     # gives histogram aspect to the bar diagram

        ax = plt.axes()
        ax.set_xticks(pos + (width / 2))
        ax.set_xticklabels(IUPAC.protein.letters)

        plt.bar(pos,
                [frequencies[letter] for letter in IUPAC.protein.letters],
                width, color='r')
        frequencies_path = os.path.join(
            freq_dir_path, "frequencies_{0}.png".format(cluster_id))
        fig.savefig(frequencies_path)


def gap_free_hamming(s, t):
    return sum(int(s[i] != t[i] and s[i] != '-' and t[i] != '-') for i in range(len(s)))


def reformat_chain(chains):
    seqs = [chains.get(key) for key in chains.keys()]
    keys = {key: i for i, key in enumerate(chains.keys())}
    return keys, seqs


def reformat_chain_dict((keys, seqs)):
    return {key: seqs[keys[key]] for key in keys}


def make_duplicates((keys, seqs), prct):
    if not prct:
        return keys, seqs
    logging.debug("Removing duplicates")
    for i, s in enumerate(seqs):
        for j, t in enumerate(seqs):
            if s.count('-') > t.count('-'):
                continue
            if gap_free_hamming(s, t) < (1 - prct) * len(s):
                for key in keys:
                    if keys[key] == j:
                        keys[key] = i
            logging.debug(str(keys))
            logging.debug(str(len(seqs)))
    return keys, seqs


def make_document(clustername, fasta):
    template = """
<!DOCTYPE html>
<html>
<head>
    <title>%s</title>
</head>
<style>
%s
</style>
<body>
<table>
%s
</table>
</body>
</html>
    """

    lines = []
    for key in fasta.keys():
        lines.append("<tr><td><small><pre>{0}</pre></small></td>"
                     "<td><small><pre>{1}</pre></small></td></tr>".format(key, color_seq(fasta.get(key))))
    return template % (clustername, style, "\n".join(lines))


def make_alignments(clusterdir):
    for cluster in os.listdir(clusterdir):
        logging.debug("cluster: {}".format(cluster))
        fasta = read_fasta(os.path.join(clusterdir, cluster), False)
        pos = cluster.rfind('.')
        name = cluster[:pos]
        with open(os.path.join(clusterdir, "{}.htm".format(name)), "wt") as fd:
            fd.write(make_document(name, fasta))


def process_one(analysis_dir, chain, chains, prct=False):
    chains = reformat_chain_dict(make_duplicates(reformat_chain(chains), prct))
    #chains = reformat_chain_dict(reformat_chain(chains))
    write_fasta(os.path.join(analysis_dir, "{}.fa".format(chain)), chains)
    processing(os.path.join(analysis_dir, "{}.fa".format(chain)),
               os.path.join(analysis_dir, chain), prct)
    make_alignments(os.path.join(os.path.join(analysis_dir, chain), "clusters"))


def process(directory, analysis_name="analysis", prct=False):
    analysis_dir = os.path.join(directory, analysis_name)
    if os.path.isdir(analysis_dir):
        shutil.rmtree(analysis_dir)
    os.mkdir(analysis_dir)

    common = get_common_name(directory, "-amino-aligned.fa")
    process_one(analysis_dir, "light", read_fasta(os.path.join(directory, common + "light-chains-amino-aligned.fa"),
                                                  False), prct)
    process_one(analysis_dir, "heavy", read_fasta(os.path.join(directory, common + "heavy-chains-amino-aligned.fa"),
                                                  False), prct)