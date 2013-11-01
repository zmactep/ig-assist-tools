__author__ = 'yakovlev'

from collections import Counter

style = """
.A, .V, .F, .P, .M, .I, .L, .W {
    color: red !important;
}

.D, .E {
    color: blue !important;
}

.R, .K {
    color: fuchsia !important;
}

.S, .T, .Y, .H, .C, .N, .G, .Q {
    color: green !important;
}

pre {
    font-family: monospace, monospace;
    font-family: 'courier new', monospace;
    font-size: 6;
}
"""


def color_seq(seq):
    return "".join("<span class=\"{0}\">{1}</span>".format(s, s) for s in seq)


def get_intersect(hcl, lcl):
    result = []
    for h in hcl.keys():
        if h[:-3]+"-VL" in lcl:
            result.append(h[:-3])
    return result


def filter_clone_name(clone_name):
    cl = clone_name[clone_name.rfind("_"):]
    if cl.startswith("_amino") or cl.startswith("_nucleo"):
        return clone_name[:clone_name.rfind("_")]
    return clone_name