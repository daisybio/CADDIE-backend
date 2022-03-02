import graph_tool as gt
import graph_tool.topology as gtt
import itertools as it


def steiner_tree(g, seeds, seed_map, weights, non_zero_hub_penalty):
    """
    :param g: graph tool Graph object
    :param seeds: list of node names of g
    :param seed_map: {graph_id : graph index}
    :param weights:
    :param non_zero_hub_penalty:
    :return:
    """

    # mc??
    mc = gt.Graph(directed=False)
    eprop_dist = mc.new_edge_property("int")
    mc.ep['dist'] = eprop_dist
    vprop_name = mc.new_vertex_property("string")
    mc.vp['graphId'] = vprop_name

    eprop_path = mc.new_edge_property("object")
    mc.ep['path'] = eprop_path
    # {gene name: vertex object}
    mc_vertex_map = dict()
    # {vertex object: seeds index}
    mc_id_map = dict()
    for i in range(len(seeds)):
        vert = mc.add_vertex()
        vprop_name[i] = seeds[i]
        mc_vertex_map[seeds[i]] = vert
        mc_id_map[vert] = i

    # the following step takes very long
    for u, v in it.combinations(seeds, 2):
        _, elist = gtt.shortest_path(g, g.vertex(seed_map[u]), g.vertex(seed_map[v]), weights=weights, negative_weights=False, pred_map=None, dag=False)
        e = mc.add_edge(mc_vertex_map[u], mc_vertex_map[v])
        eprop_dist[e] = len(elist)
        mc.ep.path[e] = list(elist)

    mst = gtt.min_spanning_tree(mc, weights=eprop_dist, root=None, tree_map=None)
    mc.set_edge_filter(mst)

    g2 = gt.Graph(directed=False)
    vprop_name = g2.new_vertex_property("string")
    g2.vp['graphId'] = vprop_name

    g2_vertex_map = dict()
    g2_id_map = dict()
    addedNodes = set()
    for i in range(len(seeds)):
        vert = g2.add_vertex()
        vprop_name[i] = seeds[i]
        g2_vertex_map[seeds[i]] = vert
        g2_id_map[vert] = i
        addedNodes.add(seeds[i])

    allmcedges = []
    for mc_edges in mc.edges():
        path = mc.ep.path[mc_edges]
        allmcedges.extend(path)
    j = len(seeds)

    allmcedges_g2 = []
    for e in allmcedges:
        sourceName = g.vertex_properties["graphId"][e.source()]
        targetName = g.vertex_properties["graphId"][e.target()]
        if sourceName not in addedNodes:
            vert = g2.add_vertex()
            vprop_name[j] = sourceName
            g2_vertex_map[sourceName] = vert
            g2_id_map[vert] = j
            addedNodes.add(sourceName)
            j += 1
        if targetName not in addedNodes:
            vert = g2.add_vertex()
            vprop_name[j] = targetName
            g2_vertex_map[targetName] = vert
            g2_id_map[vert] = j
            addedNodes.add(targetName)
            j += 1
        allmcedges_g2.append(g2.add_edge(g2_vertex_map[sourceName], g2_vertex_map[targetName]))
    weights_g2 = g2.new_edge_property("double", val=1.0)
    if non_zero_hub_penalty:
        for e, e_g2 in zip(allmcedges, allmcedges_g2):
            weights_g2[e_g2] = weights[e]
    mst2 = gtt.min_spanning_tree(g2, root=None, tree_map=None, weights=weights_g2)
    g2.set_edge_filter(mst2)
    # vw = gt.GraphView(g2, efilt=mst2)
    # g3 = Graph(vw, prune=True)

    # g3 = Graph(g22)
    while True:
        noneSteinerLeaves = []
        for i in range(g2.num_vertices()):
            if g2.vertex(i).out_degree() == 1 and g2.vertex_properties["graphId"][i] not in seeds:
                noneSteinerLeaves.append(i)
        if len(noneSteinerLeaves) == 0:
            break
        noneSteinerLeaves = reversed(sorted(noneSteinerLeaves))
        for node in noneSteinerLeaves:
            # outarray = g3.get_out_edges(node)
            g2.remove_edge(g2.edge(g2.vertex(node), g2.get_all_neighbors(node)[0]))
            g2.remove_vertex(node)
            
    return g2



