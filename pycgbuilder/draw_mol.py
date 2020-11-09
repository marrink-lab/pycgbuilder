
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import math
import networkx as nx


def rot(x, y, theta):
    return x*math.cos(theta) - y*math.sin(theta), x*math.sin(theta) + y*math.cos(theta)


def make_edge(p0, p1, order, spacing=0.1, sep=0.2):
    # Step 1: Make a horizontal line from sep to bond_length - sep, and figure
    #         out how far we'll need to rotate it. This line is (x0p, x1p)
    x0, y0 = p0
    x1, y1 = p1
    bond_length = math.sqrt((x0 - x1)**2 + (y0-y1)**2)
    sep = min(sep * bond_length, sep)
    x0p = sep
    x1p = bond_length - sep
    theta = math.atan2((y1 - y0), (x1 - x0))

    # Step 2: Split the line (if needed) in multiple copies (let's say 2). 1 up,
    #         1 down. This results in horizontal line(s) ((x0p, y0p), (x1p, y1p)).
    #         These we then rotate , resulting in ((x0pp, y0pp), (x1pp, y1pp)).
    #         Finally, translate them from the origin to (x0, y0)
    out = []
    for f in range(order):
        y0p = y1p = spacing * ((1-order)/2 + f)
        x0pp, y0pp = rot(x0p, y0p, theta)
        x1pp, y1pp = rot(x1p, y1p, theta)
        x0pp += x0
        y0pp += y0
        x1pp += x0
        y1pp += y0
        out.append([[x0pp, y0pp], [x1pp, y1pp]])

    return out


def draw_molecule(graph, clusters=None, labels=None, edge_widths=None, pos=None, ax=None):
    if not ax:
        ax = plt.gca()

    if labels is None:
        elems = nx.get_node_attributes(graph, 'element')
        labels = {}
        for n_idx, elem in elems.items():
            labels[n_idx] = '${}_{{{}}}$'.format(elem, n_idx)

    if not pos:
        pos = nx.kamada_kawai_layout(graph)

    for idx, label in labels.items():
        x, y = pos[idx]
        ax.text(x, y, label, verticalalignment='center_baseline', horizontalalignment='center')

    edges = []
    arom_edges = []
    for idx, jdx, order in graph.edges(data='order'):
        if not order:
            order = 1
        if order == 1.5:
            tmp = make_edge(pos[idx], pos[jdx], 2)
            arom_edges.append(tmp[0])
            edges.append(tmp[1])
        else:
            edges.extend(make_edge(pos[idx], pos[jdx], order))
    ax.add_collection(LineCollection(edges, color='black', linewidths=edge_widths))
    ax.add_collection(LineCollection(arom_edges, color='black', linestyle='dotted', linewidths=edge_widths))

    if clusters is not None:
        cmap = plt.get_cmap('tab10')
        for idx, position in pos.items():
            member = clusters[idx]
            ax.pie(member, center=position, radius=0.45,
                    colors=[cmap.colors[idx%cmap.N] for idx in range(len(member))])
