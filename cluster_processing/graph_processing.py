__author__ = 'yakovlev'

import os
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from Bio.SeqUtils.CheckSum import seguid

from ighumanizer3.extra.share.fasta_tools import read_fasta


def domain_name(name):
    return name[:-3]


def load_clusters(directory, strtype):
    result = []
    clusters_dir = os.path.join(os.path.join(directory, strtype), "clusters")
    for cluster in os.listdir(clusters_dir):
        result.append(read_fasta(os.path.join(clusters_dir, cluster)))
    return result


def in_one_cluster(list1, list2, cluster_list):
    for i1 in list1:
        for i2 in list2:
            for cluster in cluster_list:
                if i1 in cluster.keys() and i2 in cluster.keys():
                    return True
                break
    return False


def get_node_lists(chains):
    result = {}
    for name in chains.keys():
        seq = seguid(chains.get(name))
        if seq in result:
            result[seq].append(name)
        else:
            result[seq] = [name]
    return list(result.values())


def color_nodes(node_list, cluster_list):
    ll = len(node_list)
    #random_colors = list(np.random.random((ll, 3)))
    random_colors = list(range(ll))
    for i in node_list:
        for j in node_list:
            if in_one_cluster(i, j, cluster_list):
                random_colors[j] = random_colors[i]
    return random_colors


def is_connected(hnode, lnode):
    for h in hnode:
        for l in lnode:
            if domain_name(h) == domain_name(l):
                return True
    return False


def get_edges_list(heavy_nodes, light_nodes):
    edges = []
    for i, h in enumerate(heavy_nodes):
        for j, l in enumerate(light_nodes):
            if is_connected(h, l):
                edges.append((i, j))
    return edges


def get_labels(nodes):
    return ["\n".join(lab) for lab in nodes]


def unify_pair(heavy, light):
    result = [i for i in heavy]
    for i in light:
        result.append(i)
    return result


def unify_edges(shift, edges, labels):
    result = []
    for i, j in edges:
        result.append((labels[i], labels[shift + j]))
    return result


def get_vertices(edges):
    result = []
    for i, j in edges:
        if i not in result:
            result.append(i)
        if j not in result:
            result.append(j)
    return result


def draw_light_heavy_graph(directory):
    analysis_path = os.path.join(directory, "analysis")

    heavy_path = os.path.join(analysis_path, "heavy.fa")
    light_path = os.path.join(analysis_path, "light.fa")

    heavy = read_fasta(heavy_path)
    light = read_fasta(light_path)

    heavy_clusters = load_clusters(analysis_path, "heavy")
    light_clusters = load_clusters(analysis_path, "light")

    heavy_nodes = get_node_lists(heavy)
    light_nodes = get_node_lists(light)

    heavy_colors = color_nodes(heavy_nodes, heavy_clusters)
    light_colors = color_nodes(light_nodes, light_clusters)
    colors = unify_pair(heavy_colors, light_colors)

    heavy_labels = get_labels(heavy_nodes)
    light_labels = get_labels(light_nodes)
    labels = unify_pair(heavy_labels, light_labels)

    hl_edges = get_edges_list(heavy_nodes, light_nodes)
    edges = unify_edges(len(heavy_nodes), hl_edges, labels)

    G = nx.Graph()
    G.add_nodes_from(get_vertices(edges))
    G.add_edges_from(edges)
    nx.draw_networkx(G, color=colors)
    plt.savefig(os.path.join(analysis_path, "connections.pdf"))