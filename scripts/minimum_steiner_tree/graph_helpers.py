import numpy as np
from collections import defaultdict

from graph_tool import GraphView
from graph_tool.search import BFSVisitor


def extract_edges_from_pred(source, target, pred):
    """edges from `target` to `source` using predecessor map, `pred`"""
    edges = []
    c = target
    while c != source and pred[c] != -1:
        edges.append((pred[c], c))
        c = pred[c]
    return edges


class DistPredVisitor(BFSVisitor):
    """visitor to track distance and predecessor"""

    def __init__(self, pred, dist):
        """np.ndarray"""
        self.pred = pred
        self.dist = dist

    def tree_edge(self, e):
        s, t = int(e.source()), int(e.target())
        self.pred[t] = s
        self.dist[t] = self.dist[s] + 1


def init_visitor(g, root):
    dist = defaultdict(lambda: -1)
    dist[root] = 0
    pred = defaultdict(lambda: -1)
    vis = DistPredVisitor(pred, dist)
    return vis


def is_tree(t):
    # to undirected
    t = GraphView(t, directed=False)

    # num nodes = num edges+1
    if t.num_vertices() != (t.num_edges() + 1):
        return False

    # all nodes have degree > 0
    vs = list(map(int, t.vertices()))
    degs = t.degree_property_map('out').a[vs]
    if np.all(degs > 0) == 0:
        return False

    return True


def is_steiner_tree(t, X):
    if not is_tree(t):
        return False
    for x in X:
        if not has_vertex(t, x):
            return False
    return True


def has_vertex(g, i):
    # to avoid calling g.vertex
    return g._Graph__filter_state['vertex_filter'][0].a[i] > 0
