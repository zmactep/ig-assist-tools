__author__ = 'yakovlev'

import os
import shutil

from ighumanizer3.extra.share.fasta_tools import read_fasta, fasta_from_keylist, write_fasta
from ighumanizer3.extra.share.algorithm import get_common_name
from similarity_processing.report_generator import generate_report


def cluster_projname(cluster):
    pattern = cluster[0].split('_')[0]
    if all(name.split('_')[0] == pattern for name in cluster):
        return pattern
    else:
        return "Merged"


def is_similar(s, t, shead, prct):
    ss = s.rstrip('-')
    tt = t.rstrip('-')
    lL = max((len(ss), ss), (len(tt), tt))
    ll = min((len(ss), ss), (len(tt), tt))
    lL, ll = lL[1], ll[1]
    if len(ll) < len(lL) * prct:
        return False
    # not only gaps
    # for i in range(ll[0]):
    #     if s[i] != t[i]:
    #         return False
    # return True
    return lL[shead:].startswith(ll[shead:int(len(lL) * prct)])


def get_sclusters(chains, minlen, shead, prct):
    clusters = []
    others = []
    visited = {key: False for key in chains.keys()}
    for key in chains.keys():
        # Already in cluster
        if visited[key]:
            continue
        # Bad chain
        if len(chains[key].replace('-', '')) < minlen:
            others.append(key)
            visited[key] = True
        # Make new cluster
        visited[key] = True
        clusters.append([key])
        # Check all unvisited
        for other_key in filter(lambda k: not visited[k], visited.keys()):
            # Bad chains
            if len(chains[other_key].replace('-', '')) < minlen:
                others.append(other_key)
                visited[other_key] = True
            # Append to cluster
            if is_similar(chains[key], chains[other_key], shead, prct):
                clusters[-1].append(other_key)
                visited[other_key] = True
    clusters_d = {"%s_G%i" % (cluster_projname(c), i): fasta_from_keylist(chains, c) for i, c in enumerate(clusters)}
    if others:
        clusters_d["%s_others" % cluster_projname(others)] = fasta_from_keylist(chains, others)
    return clusters_d


def process_one(similarity_dir, chain, chains, minlen,  shead, prct):
    writepath = os.path.join(similarity_dir, chain)
    os.mkdir(writepath)
    clusters = get_sclusters(chains, minlen, shead, prct)
    for key in clusters:
        write_fasta(os.path.join(writepath, key + ".fa"), clusters[key])
    group_fasta = {k: max((clusters[k][name] for name in clusters[k]), key=len) for k in clusters.keys()}
    write_fasta(os.path.join(similarity_dir, "%s-groups.fa" % chain), group_fasta)
    return clusters


def process(directory, similarity_name="similarity", minlen=100, shead=0, prct=1.):
    similarity_dir = os.path.join(directory, similarity_name)
    if os.path.isdir(similarity_dir):
        shutil.rmtree(similarity_dir)
    os.mkdir(similarity_dir)

    common = get_common_name(directory, "-amino-aligned.fa")
    print(common)
    light = process_one(similarity_dir, "light", read_fasta(os.path.join(directory,
                                                                         common + "light-chains-amino-aligned.fa"),
                                                            False), minlen, shead, prct)
    heavy = process_one(similarity_dir, "heavy", read_fasta(os.path.join(directory,
                                                                         common + "heavy-chains-amino-aligned.fa"),
                                                            False), minlen, shead, prct)

    with open(os.path.join(similarity_dir, "report.htm"), "wt") as fd:
        fd.write(generate_report(heavy, light, minlen, shead, prct))