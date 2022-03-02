import sys

__time = 0


def __dfs_find_bridges(g, node, visited, disc, low, parent, is_bridge):
    
    visited[node] = True
    global __time
    disc[node] = __time
    low[node] = __time
    __time += 1
    
    for nb in g.get_all_neighbors(node):
        if not visited[nb]:
            parent[nb] = node
            __dfs_find_bridges(g, int(nb), visited, disc, low, parent, is_bridge)
            low[node] = min(low[node], low[nb])
            if low[nb] > disc[node]:
                is_bridge[g.edge(node, nb)] = True
        elif int(nb) != parent[node]:
            low[node] = min(low[node], disc[nb])

def find_bridges(g):
    r"""Finds all bridges in a graph."""
    
    global __time
    __time = 0
    sys.setrecursionlimit(g.num_vertices() + 1)
    visited = g.new_vertex_property("boolean", False)
    disc = g.new_vertex_property("float", float("inf"))
    low = g.new_vertex_property("float", float("inf"))
    parent = g.new_vertex_property("int", -1)
    is_bridge = g.new_edge_property("boolean", False)
    for node in range(g.num_vertices()):
        if not visited[node]:
            __dfs_find_bridges(g, node, visited, disc, low, parent, is_bridge)
    return is_bridge
