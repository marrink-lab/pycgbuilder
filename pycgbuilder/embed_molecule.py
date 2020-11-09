import networkx as nx
import numpy as np
import math


def rescale_bondlengths(graph, embedding, scale=1, reductor=np.median):
    keys = list(embedding)
    positions = [embedding[key] for key in keys]
    positions = np.array(positions)
    bonds = []
    for idx, jdx in graph.edges:
        length = np.linalg.norm(positions[keys.index(idx)] - positions[keys.index(jdx)])
        bonds.append(length)
    min_len = reductor(bonds)
    mean = np.mean(positions, axis=0)
    positions -= mean
    positions *= scale/min_len
    positions += mean
    return {key: positions[idx] for idx, key in enumerate(keys)}


def rescale(layout_func):
    def layout(graph, *args, **kwargs):
        try:
            pos = layout_func(graph, *args, **kwargs)
        except:
            pos = {}
        return rescale_bondlengths(graph, pos)
    return layout


@rescale
def vsepr_layout(graph):
    dist = dict(nx.shortest_path_length(graph, weight=None))
    for source, dest_d in dist.items():
        for dest, d in dest_d.items():
            if d % 2 == 1:  # odd
                # sqrt(1 + (n*x)**2 + n*x**2)
                n = (d + 1)/2
                dist[source][dest] = math.sqrt(3*n**2 - 3*n + 1)
            else:  # even
                dist[source][dest] = d / 2 * math.sqrt(3)
    try:
        p0 = nx.planar_layout(graph)
    except nx.NetworkXException:
        p0 = None
    pos = nx.kamada_kawai_layout(graph, pos=p0, weight=None, dist=dist, scale=len(graph) ** 0.5)
    return pos


kamada_kawai_layout = rescale(nx.kamada_kawai_layout)
spring_layout = rescale(nx.spring_layout)
spectral_layout = rescale(nx.spectral_layout)
planar_layout = rescale(nx.planar_layout)
