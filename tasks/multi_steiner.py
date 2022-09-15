from tasks.task_hook import TaskHook
from tasks.util.steiner_tree import steiner_tree
from tasks.util.find_bridges import find_bridges
from tasks.util.read_graph_tool_graph import read_graph_tool_graph
from tasks.util.edge_weights import edge_weights
import os.path
import graph_tool as gt
import graph_tool.util as gtu
import graph_tool.topology as gtt
import sys
from caddie import models
import copy


def multi_steiner(task_hook: TaskHook):

    seeds = task_hook.parameters["seeds"]
    seeds.sort()

    target = task_hook.parameters.get("target", "drug-target")

    cancer_types = task_hook.parameters.get("cancer_types", [])
   
    cancer_dataset = task_hook.parameters.get("cancer_dataset", "NCG6")

    gene_interaction_datasets = task_hook.parameters.get("gene_interaction_datasets", ["BioGRID"])
    
    drug_interaction_datasets = task_hook.parameters.get("drug_interaction_datasets", ["BioGRID"])
    
    ignored_edge_types = task_hook.parameters.get("ignored_edge_types", [])
    
    ignore_non_seed_baits = task_hook.parameters.get("ignore_non_seed_baits", False)

    num_trees = task_hook.parameters.get("num_trees", 5)
    
    tolerance = task_hook.parameters.get("tolerance", 10)
    
    max_deg = task_hook.parameters.get("max_deg", sys.maxsize)
    
    hub_penalty = task_hook.parameters.get("hub_penalty", 0.0)
    
    num_threads = task_hook.parameters.get("num_threads", 1)

    mutation_cancer_type = task_hook.parameters.get("mutation_cancer_type", None)
    if mutation_cancer_type is not None:
      mutation_cancer_type = models.MutationCancerType.objects.filter(name__iexact=mutation_cancer_type).first()
      if mutation_cancer_type is None:
        raise ValueError('Could not find mutation_cancer_type.')
      
    expression_cancer_type = task_hook.parameters.get("expression_cancer_type", None)
    if expression_cancer_type is not None:
      expression_cancer_type = models.ExpressionCancerType.objects.filter(name__iexact=expression_cancer_type).first()
      if expression_cancer_type is None:
        raise ValueError('Could not find expression_cancer_type.')
    
    tissue = task_hook.parameters.get("tissue", None)
    if tissue is not None:
      tissue = models.Tissue.objects.filter(name__iexact=tissue).first()
      if tissue is None:
        raise ValueError('Could not find tissue.')

    # Set number of threads if OpenMP support is enabled.
    if gt.openmp_enabled():
        gt.openmp_set_num_threads(num_threads)
    
    # Parsing input file.
    task_hook.set_progress(0 / (float(num_trees + 3)), "Parsing input.")
    file_path = os.path.join(
        task_hook.data_directory,
        f"internal_{cancer_dataset}.gt"    )

    # seeds ids is list of graph ids of seed nodes
    g, seed_ids, _, drug_ids, degrees = read_graph_tool_graph(
        file_path=file_path,
        gene_datasets=gene_interaction_datasets,
        drug_datasets=drug_interaction_datasets,
        seeds=seeds,
        cancer_types=cancer_types,
        ignored_edge_types=ignored_edge_types,
        max_deg=max_deg,
        ignore_non_seed_baits=ignore_non_seed_baits,
        target=target
    )

    # map is {node graphId : network id}
    seed_map = {g.vertex_properties['graphId'][node]: node for node in seed_ids}

    task_hook.set_progress(1 / (float(num_trees + 3)), "Computing edge weights.")
    weights = edge_weights(
        g,
        hub_penalty,
        mutation_cancer_type,
        expression_cancer_type,
        tissue,
        inverse=True,
    )

    # Find first steiner trees
    task_hook.set_progress(2 / (float(num_trees + 3)), "Computing Steiner tree 1 of {}.".format(num_trees))
    first_tree = steiner_tree(g, seeds, seed_map, weights, hub_penalty > 0)

    num_found_trees = 1
    tree_edges = []

    for tree_edge in first_tree.edges():
        source_name = first_tree.vertex_properties["graphId"][first_tree.vertex_index[tree_edge.source()]]
        target_name = first_tree.vertex_properties["graphId"][first_tree.vertex_index[tree_edge.target()]]
        tree_edges.append((gtu.find_vertex(g, prop=g.vertex_properties['graphId'], match=source_name)[0], gtu.find_vertex(g, prop=g.vertex_properties['graphId'], match=target_name)[0]))
    cost_first_tree = sum([weights[g.edge(source, target)] for source, target in tree_edges])

    returned_nodes = set(int(gtu.find_vertex(g, prop=g.vertex_properties['graphId'], match=first_tree.vertex_properties['graphId'][node])[0]) for node in range(first_tree.num_vertices()))
    if num_trees > 1:
        is_bridge = find_bridges(g)
        edge_filter = g.new_edge_property("boolean", True)
        found_new_tree = True
        while len(tree_edges) > 0:
            if found_new_tree:
                task_hook.set_progress(float(num_found_trees + 2) / (float(num_trees + 3)), "Computing Steiner tree {} of {}.".format(num_found_trees + 1, num_trees))
            found_new_tree = False
            tree_edge = tree_edges.pop()
            g_edge = g.edge(tree_edge[0], tree_edge[1])
            if not is_bridge[g_edge]:
                edge_filter[g_edge] = False
                g.set_edge_filter(edge_filter)
                next_tree = steiner_tree(g, seeds, seed_map, weights, hub_penalty > 0)
                next_tree_edges = set()
                for next_tree_edge in next_tree.edges():
                    source_name = next_tree.vertex_properties["graphId"][next_tree.vertex_index[next_tree_edge.source()]]
                    target_name = next_tree.vertex_properties["graphId"][next_tree.vertex_index[next_tree_edge.target()]]
                    next_tree_edges.add((gtu.find_vertex(g, prop=g.vertex_properties['graphId'], match=source_name)[0], gtu.find_vertex(g, prop=g.vertex_properties['graphId'], match=target_name)[0]))
                cost_next_tree = sum([weights[g.edge(source, target)] for source, target in next_tree_edges])
                if cost_next_tree <= cost_first_tree * ((100.0 + tolerance) / 100.0):
                    found_new_tree = True
                    num_found_trees += 1
                    for node in range(next_tree.num_vertices()):
                        returned_nodes.add(int(gtu.find_vertex(g, prop=g.vertex_properties['graphId'], match=next_tree.vertex_properties["graphId"][node])[0]))
                    removed_edges = []
                    for source, target in tree_edges:
                        if not ((source, target) in set(next_tree_edges)) or ((target, source) in set(next_tree_edges)):
                            removed_edges.append((source, target))
                    for edge in removed_edges:
                        tree_edges.remove(edge)
                g.clear_filters()
                edge_filter[g_edge] = True
            if num_found_trees >= num_trees:
                break

    task_hook.set_progress((float(num_trees + 2)) / (float(num_trees + 3)), "Formatting results")

    returned_edges = []
    returned_nodes_with_connectors = copy.deepcopy(seed_ids) # return seed_ids in any case, even if no edges
    for node_i in returned_nodes:
        for node_j in returned_nodes:
            if node_i == node_j:
                continue
            # not from drug to drug
            if g.vertex_properties["type"][node_i] == g.vertex_properties["type"][node_j] == 'Drug':
                continue

            vertices, edges = gtt.shortest_path(g, node_i, node_j)

            # if drug node is in path, skip
            # we dont want to add paths over drugs
            # Should have no effect since drugs are filtered out of graph? TODO
            continue_flag = False
            for vertex in vertices:
                if g.vertex_properties["type"][int(vertex)] == 'Drug' and int(vertex) not in returned_nodes:
                    continue_flag = True
            if continue_flag:
                continue

            for vertex in vertices:
                if int(vertex) not in returned_nodes_with_connectors:
                    returned_nodes_with_connectors.append(int(vertex))
            for edge in edges:
                if ((edge.source(), edge.target()) not in returned_edges) or ((edge.target(), edge.source()) not in returned_edges):
                    returned_edges.append((edge.source(), edge.target()))

    subgraph = {"nodes": [g.vertex_properties["graphId"][node] for node in returned_nodes_with_connectors],
                "edges": [{"from": g.vertex_properties["graphId"][source], "to": g.vertex_properties["graphId"][target]} for source, target in returned_edges]}
    node_types = {g.vertex_properties["graphId"][node]: g.vertex_properties["type"][node] for node in returned_nodes_with_connectors}
    is_seed = {g.vertex_properties["graphId"][node]: node in set(seed_ids) for node in returned_nodes_with_connectors}
    is_result = {g.vertex_properties["graphId"][node]: node in returned_nodes for node in returned_nodes_with_connectors}
    db_degrees = {g.vertex_properties["graphId"][node]: degrees[g.vertex_properties["graphId"][node]] for node in returned_nodes}

    task_hook.set_results({
        "network": subgraph,
        "node_attributes": {"node_types": node_types, "is_seed": is_seed, "is_result": is_result, "db_degrees": db_degrees},
        "cancer_types": cancer_types,
        'cancer_dataset': cancer_dataset,
        'gene_interaction_dataset': gene_interaction_datasets[0],
        'drug_interaction_dataset': drug_interaction_datasets[0]
    })
