import networkx as nx
import math


def vsepr_layout(graph):
    dist = dict(nx.shortest_path_length(graph, weight=None))
    for source, dest_d in dist.items():
        for dest, d in dest_d.items():
            if d % 2 == 1:  # odd
                # sqrt(1 + (n*x)**2 + n*x**2)
                n = (d - 1) // 2
                dist[source][dest] = math.sqrt(1 + 3 * n ** 2 + 3 * n)
            else:  # even
                dist[source][dest] = d // 2 * math.sqrt(3)
    try:
        p0 = nx.planar_layout(graph)
    except nx.NetworkXException:
        p0 = None
    pos = nx.kamada_kawai_layout(graph, pos=p0, weight=None, dist=dist, scale=len(graph) ** 0.5)
    return pos


# vsepr_layout = nx.kamada_kawai_layout
