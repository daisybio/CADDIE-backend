import graph_tool.topology as gtt

def scores_to_results(
        target,
        result_size,
        g,
        seed_ids,
        cancer_driver_gene_ids,
        drug_ids,
        scores,
        cancer_types,
        cancer_dataset,
        gene_interaction_dataset,
        drug_interaction_dataset,
        filterPaths,
        degrees
):
    r"""Transforms the scores to the required result format."""
    candidates = []
    if target == "drug":
        candidates = [(node, scores[node]) for node in drug_ids if scores[node] > 0]
    else:
        candidates = [(node, scores[node]) for node in range(g.num_vertices()) if scores[node] > 0 and node not in set(seed_ids)]
    best_candidates = [item[0] for item in sorted(candidates, key=lambda item: item[1], reverse=True)[:result_size]]
    # Concatenate best result candidates with seeds and compute induced subgraph.
    # since the result size filters out nodes, we have to find them here again to compute the complete result network
    # in case there is no direct interaction

    returned_edges = set()
    returned_nodes = set(seed_ids) # return seed_ids in any case

    # return only the path to a drug with the shortest distance
    if filterPaths:
        for candidate in best_candidates:
            distances = gtt.shortest_distance(g, candidate, seed_ids)
            closest_distance_mean = sum(distances) / len(distances)

            for index, seed_id in enumerate(seed_ids):
                if distances[index] > closest_distance_mean:
                    continue
                vertices, edges = gtt.shortest_path(g, candidate, seed_id)

                drug_in_path = False
                for vertex in vertices:
                    if g.vertex_properties["type"][int(vertex)] == "Drug" and vertex != candidate:
                        drug_in_path = True
                        break
                if drug_in_path:
                    continue

                for vertex in vertices:
                    if int(vertex) not in returned_nodes:
                        returned_nodes.add(int(vertex))
                for edge in edges:
                    if ((edge.source(), edge.target()) not in returned_edges) or ((edge.target(), edge.source()) not in returned_edges):
                        returned_edges.add((edge.source(), edge.target()))
    else:
        for candidate in best_candidates:
            for index, seed_id in enumerate(seed_ids):
                vertices, edges = gtt.shortest_path(g, candidate, seed_id)

                drug_in_path = False
                for vertex in vertices:
                    if g.vertex_properties["type"][int(vertex)] == "Drug" and vertex != candidate:
                        drug_in_path = True
                        break
                if drug_in_path:
                    continue

                for vertex in vertices:
                    if int(vertex) not in returned_nodes:
                        returned_nodes.add(int(vertex))
                for edge in edges:
                    if ((edge.source(), edge.target()) not in returned_edges) or ((edge.target(), edge.source()) not in returned_edges):
                        returned_edges.add((edge.source(), edge.target()))


            # for neighbor in g.get_all_neighbors(node):
            #     if int(neighbor) > node and int(neighbor) in returned_nodes:
            #         returned_edges.append((node, int(neighbor)))
    # max_score = max([scores[node] for node in returned_nodes])
    # max_non_virus_score = max([scores[node] for node in returned_nodes if not node in set(cancer_driver_gene_ids)])
    # max_non_seed_score = float("nan")
    # if result_size > 0:
    #     max_non_seed_score = max([scores[node] for node in best_candidates])
    subgraph = {"nodes": [g.vertex_properties["graphId"][node] for node in returned_nodes],
                "edges": [{"from": g.vertex_properties["graphId"][source], "to": g.vertex_properties["graphId"][target]} for source, target in returned_edges],
                # "max_score": max_score,
                # "max_non_virus_score": max_non_virus_score,
                # "max_non_seed_score": max_non_seed_score,
                }

    # Compute node attributes.
    node_types = {g.vertex_properties["graphId"][node]: g.vertex_properties["type"][node] for node in returned_nodes}
    is_seed = {g.vertex_properties["graphId"][node]: node in set(seed_ids) for node in returned_nodes}
    returned_scores = {g.vertex_properties["graphId"][node]: scores[node] for node in returned_nodes}
    is_result = {g.vertex_properties["graphId"][node]: node in best_candidates for node in returned_nodes}
    db_degrees = {g.vertex_properties["graphId"][node]: degrees[g.vertex_properties["graphId"][node]] for node in returned_nodes}

    traces_degree = {'Drug': {'x': [], 'y': [], 'names': []}, 'CancerNode': {'x': [], 'y': [], 'names': []}, 'Node': {'x': [], 'y': [], 'names': []}}
    for node in returned_nodes:
        degree = degrees[g.vertex_properties["graphId"][node]]
        score = returned_scores[g.vertex_properties["graphId"][node]]
        node_type = node_types[g.vertex_properties["graphId"][node]]
        traces_degree[node_type]['x'].append(degree)
        traces_degree[node_type]['y'].append(score)
        traces_degree[node_type]['names'].append(g.vertex_properties["name"][node])

    return {
        "network": subgraph,
        "node_attributes":
            {
                "node_types": node_types,
                "is_seed": is_seed,
                "scores": returned_scores,
                "is_result": is_result,
                'db_degrees': db_degrees,
            },
        'traces_degree': traces_degree,
        'cancer_types': cancer_types,
        'cancer_dataset': cancer_dataset,
        'gene_interaction_dataset': gene_interaction_dataset,
        'drug_interaction_dataset': drug_interaction_dataset
    }
