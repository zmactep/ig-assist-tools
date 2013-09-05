__author__ = 'mactep'


def gap_free_hamming(s, t):
    return sum(int(s[i] == t[i] or s[i] == '-' or t[i] == 'i') for i in range(len(s)))


def change_key_val(keys, old, new):
    for key in keys:
        if keys[key] == old:
            keys[key] = new
    return keys


def reformat_chain(chains):
    seqs = [chains.get(key) for key in chains.keys()]
    keys = {key: i for i, key in enumerate(chains.keys())}
    return keys, seqs


def reformat_chain(keys, seqs):
    return {key: seqs[keys[key]] for key in keys}


def remove_duplicates(keys, seqs):
    for i, s in enumerate(seqs):
        for j, t in enumerate(seqs):
            if s.count('-') > t.count('-'):
                continue
            if gap_free_hamming(s, t) < 4:
                keys = change_key_val(keys, j, i)
    return keys, seqs