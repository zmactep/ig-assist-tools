__author__ = 'yakovlev'

import os
from collections import Counter
from ighumanizer3.extra.share.fasta_tools import read_fasta, get_consensus_amino

from cluster_processing import *


def load_clusters(directory, strtype):
    result = []
    clusters_dir = os.path.join(os.path.join(directory, strtype), "clusters")
    for cluster in filter(lambda x: x.endswith(".fasta") or x.endswith(".fa"), os.listdir(clusters_dir)):
        result.append(read_fasta(os.path.join(clusters_dir, cluster), False))
    return result


def get_cells(heavy, light):
    return [[map(filter_clone_name, get_intersect(hcl, lcl)) for hcl in heavy] for lcl in light]


def get_consensus(cl, prct):
    result = []
    ll = len(cl[cl.keys()[0]])
    size = len(cl.keys())
    for i in range(ll):
        letter, count = Counter(cl.get(key)[i] for key in cl.keys()).most_common(1)[0]
        if count > size * prct:
            result.append(letter)
        else:
            result.append('X')
    return ''.join(result)


def hl_process(directory, prct):
    light = load_clusters(directory, "light")
    heavy = load_clusters(directory, "heavy")
    cells = get_cells(heavy, light)

    # header
    head = ["<tr><td></td>"]
    for i in range(len(heavy)):
        head.append("<td>VH%i</td>" % (i + 1))
    head.append("</tr>")

    lines = ["".join(head)]
    for i in range(len(light)):
        line = ["<tr><td>VL%i</td>" % (i + 1)]
        for j in range(len(heavy)):
            line.append("<td><small>%s</small></td>" % ("<br>".join(cells[i][j])))
        line.append("</tr>")
        lines.append("".join(line))

    hcons = "<br>\n".join(("VH%i : <pre>%s</pre>" % (i, color_seq(s))) for i, s in enumerate(map(get_consensus_amino, heavy)))
    lcons = "<br>\n".join(("VL%i : <pre>%s</pre>" % (i, color_seq(s))) for i, s in enumerate(map(get_consensus_amino, light)))

    template = """
<!DOCTYPE html>
<html>
<head>
    <title>Final statistics</title>
</head>
<style>
%s
</style>
<body>
<p><b>Equality confidence:</b> %0.2f</p> <br>
<table border=1>
%s
</table>
<br><br>
<p><b>Consensuses:</b></p> <br>
%s
<br><br>
%s
</body>
</html>
    """ % (style, prct or 1.0, "\n".join(lines), hcons, lcons)

    with open(os.path.join(directory, "result.htm"), "wt") as fd:
        fd.write(template)