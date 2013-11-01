__author__ = 'yakovlev'

from cluster_processing import *


def get_cells(heavy, light):
    return {(h, l): map(filter_clone_name, get_intersect(heavy[h], light[l])) for h in heavy for l in light}


# def sort_key(x):
#     if "other" in x:
#         return 10050
#     pos = x.find('_')
#     return int(x[pos + 2:])

def sorted_keys(chain):
    f_sort_key = lambda x: int(x[x.find("_") + 2:]) if "other" not in x else 100500
    return sorted(chain.keys(), key=f_sort_key)


def generate_first_table(cells, heavy, light):
    # header
    head = ["<tr><td></td>"]
    heavy_sorted = sorted_keys(heavy)
    light_sorted = sorted_keys(light)
    for hkey in heavy_sorted:
        head.append("<td>VH: %s</td>" % hkey)
    head.append("</tr>")

    lines = ["".join(head)]
    for lkey in light_sorted:
        line = ["<tr><td>VL: %s</td>" % lkey]
        for hkey in heavy_sorted:
            c = cells[hkey, lkey]
            if len(c) > 1:
                print(c)
            line.append("<td><small>%s</small></td>" % ("<br>".join(c)))
        line.append("</tr>")
        lines.append("".join(line))
    return "\n".join(lines)


def get_longest(fastadict, ckeys):
    return max(((key, value) for key, value in fastadict.items() if filter_clone_name(key) in ckeys),
               key=lambda (x, y): len(y))[0]


def get_shortest(fastadict, ckeys):
    return min(((key, value) for key, value in fastadict.items() if filter_clone_name(key) in ckeys),
               key=lambda (x, y): len(y))[0]


def generate_second_table(cells, heavy, light):
    head = "<tr>"                 \
           "<th>Group</th>"       \
           "<th>VH Group</th>"    \
           "<th>VL Group</th>"    \
           "<th>Clones</th>"      \
           "<th>VH</th>" \
           "<th>VL</th>" \
           "</tr>"
    resultlist = []
    heavy_sorted = sorted_keys(heavy)
    light_sorted = sorted_keys(light)
    lines = [head]
    group_id = 0
    for hkey in heavy_sorted:
        if "other" in hkey:
            continue
        for lkey in light_sorted:
            if "other" in lkey:
                continue
            c = cells[hkey, lkey]
            if not c:
                continue
            group_id += 1
            longest = get_longest(heavy[hkey], c)[:-3]
            resultlist.append(filter_clone_name(longest))
            line = "<tr>"                   \
                   "<td>{}</td>"            \
                   "<td>{}</td><td>{}</td>" \
                   "<td>{}</td>"            \
                   "<td><pre>{}</pre></td><td><pre>{}</pre></td>" \
                   "</tr>".format(group_id,
                                  hkey, lkey,
                                  "\n".join(c),
                                  heavy[hkey][longest + "-VH"].replace('-', ''),
                                  light[lkey][longest + "-VL"].replace('-', ''))
            lines.append(line)

    return "\n".join(lines), resultlist


def get_stats(heavy, light, cells):
    th = sum(len(heavy[i]) for i in heavy)
    tl = sum(len(light[i]) for i in light)
    notnull = sum(int(cells[i, j] != []) for i, j in cells if "others" not in i and "others" not in j)
    hothers = sum(len(heavy[i]) for i in heavy if "others" in i)
    lothers = sum(len(light[i]) for i in light if "others" in i)
    lines = ["<p><b>Total clones:</b> %i</p>" % (th if th == tl else -th),
             "<p><b>Intersections:</b> %i</p>" % notnull,
             "<p><b>Heavy others:</b> %i</p>" % hothers,
             "<p><b>Light others:</b> %i</p>" % lothers]
    return "\n".join(lines)


def generate_report(heavy, light, minlen, shead, prct):
    cells = get_cells(heavy, light)
    first_table = generate_first_table(cells, heavy, light)
    second_table, resultlist = generate_second_table(cells, heavy, light)

    template = """
<!DOCTYPE html>
<html>
<head>
    <title>Similarity statistics</title>
</head>
<style>
%s
</style>
<body>
<p><b>Minimal length:</b> %i</p>
<p><b>Head skip:</b> %i</p>
<p><b>Tail equality confidence:</b> %0.2f</p>
<br>
%s
<br>
<table border=1>
%s
</table>
<br><br>
<table border=1>
%s
</table>
<br><br>
<b>Clones:</b>
<br>
%s
</body>
</html>
    """ % (style,
           minlen, shead, prct,
           get_stats(heavy, light, cells),
           first_table, second_table,
           "<br>\n".join(resultlist))

    return template

