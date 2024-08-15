from tasks.task_hook import TaskHook
from tasks.util.read_graph_tool_graph import read_graph_tool_graph
from tasks.util.edge_weights import edge_weights, custom_edge_weights
import os.path
import graph_tool as gt
import graph_tool.topology as gtt
import sys
import numpy as np
from caddie import models


def network_proximity(task_hook: TaskHook):

    seeds = task_hook.parameters['seeds']
    
    cancer_types = task_hook.parameters.get("cancer_types", [])
    
    cancer_dataset = task_hook.parameters.get("cancer_dataset", "NCG6")
   
    gene_interaction_datasets = task_hook.parameters.get("gene_interaction_datasets", ["BioGRID"])
    
    drug_interaction_datasets = task_hook.parameters.get("drug_interaction_datasets", ["BioGRID"])
    
    include_non_approved_drugs = task_hook.parameters.get("include_non_approved_drugs", False)

    num_random_seed_sets = task_hook.parameters.get("num_random_seed_sets", 32)

    num_random_drug_target_sets = task_hook.parameters.get("num_random_drug_target_sets", 32)

    result_size = task_hook.parameters.get("result_size", 20)

    max_deg = task_hook.parameters.get("max_deg", sys.maxsize)

    hub_penalty = task_hook.parameters.get("hub_penalty", 0.0)

    num_threads = task_hook.parameters.get("num_threads", 1)

    include_nutraceutical_drugs = task_hook.parameters.get("include_nutraceutical_drugs", False)

    only_atc_l_drugs = task_hook.parameters.get("only_atc_l_drugs", False)

    include_only_ctrpv2_drugs = task_hook.parameters.get("include_only_ctrpv2_drugs", False)

    filter_paths = task_hook.parameters.get("filter_paths", True)
    
    # {'graph_id1': value1, 'graph_id2': value2, ...}
    custom_node_weights = task_hook.parameters.get("custom_node_weights", None)
    custom_node_weights_additionally = task_hook.parameters.get("custom_node_weights_additionally", False)

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

    drug_target_action = task_hook.parameters.get("drug_target_action", None)
    
    available_drugs = task_hook.parameters.get("available_drugs", None)
    if available_drugs is not None:
        available_drugs = [drug.lower() for drug in available_drugs]

    # Set number of threads if OpenMP support is enabled.
    if gt.openmp_enabled():
        gt.openmp_set_num_threads(num_threads)

    # Parsing input file.
    task_hook.set_progress(0.0 / 8, "Parsing input.")
    file_path = os.path.join(
        task_hook.data_directory,
        f"internal_{cancer_dataset}.gt"
        )
    g, seed_graph_ids, cancer_node_ids, drug_ids, degrees = read_graph_tool_graph(
        file_path=file_path,        
        gene_datasets=gene_interaction_datasets,
        drug_datasets=drug_interaction_datasets,
        seeds=seeds,
        cancer_types=cancer_types,
        ignored_edge_types=[],
        max_deg=max_deg,
        ignore_non_seed_baits=False,
        include_indirect_drugs=True,
        include_non_approved_drugs=include_non_approved_drugs,
        include_nutraceutical_drugs=include_nutraceutical_drugs,
        only_atc_l_drugs=only_atc_l_drugs,
        target='drug',
        drug_action=drug_target_action,
        available_drugs=available_drugs,
        include_only_ctrpv2_drugs=include_only_ctrpv2_drugs
    )

    # Computing edge weights.
    task_hook.set_progress(1.0 / 8, "Computing edge weights.")
    
    if custom_node_weights_additionally:
        weights = edge_weights(
            g,
            hub_penalty,
            mutation_cancer_type,
            expression_cancer_type,
            tissue,
            inverse=True,
        )
        weights = custom_edge_weights(g, custom_node_weights, weights)
    else:
        weights = custom_edge_weights(g, custom_node_weights)

    # task_hook.set_progress(2.0 / 8, "Deleting drug targets not in LCC.")
    drug_targets = {drug_id: g.get_all_neighbors(drug_id) for drug_id in drug_ids}

    # Compute all shortest path distances.
    task_hook.set_progress(3.0 / 8, "Computing all shortest path distances.")
    distances = None
    if hub_penalty == 0:
        distances = gtt.shortest_distance(g)
    else:
        distances = gtt.shortest_distance(g, weights=weights)
    # remove infinite from distances. For that we convert to numpy array
    distance_matrix = []
    for distance_vector in distances:
        vector = distance_vector.get_array()
        vector[vector == np.inf] = 9999999
        distance_matrix.append(vector)

    # Compute network proximities.
    task_hook.set_progress(4.0 / 8, "Computing network proximities.")
    proximities = {drug_id : 0 for drug_id in drug_ids}

    for drug_id in drug_ids:
        distance = 0.0
        for drug_target in drug_targets[drug_id]:
            distance += min([distance_matrix[drug_target][seed_id] for seed_id in seed_graph_ids])
        if len(drug_targets[drug_id]) == 0:
            proximities[drug_id] = np.inf
        else:
            proximities[drug_id] = distance / float(len(drug_targets[drug_id]))

    # Compute background distribution.
    task_hook.set_progress(5.0 / 8, "Computing background distribution")
    min_num_targets = min([len(drug_targets[drug_id]) for drug_id in drug_ids])
    max_num_targets = max([len(drug_targets[drug_id]) for drug_id in drug_ids])
    node_graph_ids = list(range(g.num_vertices()))
    background_distribution = []

    num_seeds = len(seed_graph_ids)
    for i in range(num_random_seed_sets):
        np.random.shuffle(node_graph_ids)
        random_seed_ids = node_graph_ids[:num_seeds]
        for k in range(num_random_drug_target_sets):
            np.random.shuffle(node_graph_ids)
            random_drug_targets = node_graph_ids[:np.random.randint(min_num_targets, max_num_targets + 1)]
            # skip calculation if 0 random drug targets in list
            if not len(random_drug_targets):
                continue
            distance = 0.0
            for drug_target in random_drug_targets:
                distance += min([distance_matrix[drug_target][seed_id] for seed_id in random_seed_ids])
            background_distribution.append(distance / float(len(random_drug_targets)))

    background_mean = np.mean(background_distribution)
    background_std = np.std(background_distribution)

    # Apply Z-score transformation.
    task_hook.set_progress(6.0 / 8, "Applying Z-score transformation.")
    drugs_with_z_scores = [(drug_id, (proximities[drug_id] - background_mean) / background_std) for drug_id in drug_ids]

    task_hook.set_progress(7.0 / 8, "Formatting results.")
    best_drugs = [item for item in sorted(drugs_with_z_scores, key=lambda item: item[1])[:result_size]]
    best_drugs_ids = [item[0] for item in best_drugs]
    seed_graph_ids = list(set(seed_graph_ids))
    returned_edges = set()

    returned_nodes = set(seed_graph_ids) # return seed_ids in any case

    if filter_paths:
        for candidate in best_drugs_ids:
            distances = gtt.shortest_distance(g, candidate, seed_graph_ids)
            closest_distance_mean = sum(distances) / len(distances)

            for index, seed_id in enumerate(seed_graph_ids):
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
                    if ((edge.source(), edge.target()) not in returned_edges) or (
                            (edge.target(), edge.source()) not in returned_edges):
                        returned_edges.add((edge.source(), edge.target()))
    else:
        for candidate in best_drugs_ids:
            for index, seed_id in enumerate(seed_graph_ids):
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
                    if ((edge.source(), edge.target()) not in returned_edges) or (
                            (edge.target(), edge.source()) not in returned_edges):
                        returned_edges.add((edge.source(), edge.target()))




    subgraph = {"nodes": [g.vertex_properties["graphId"][node] for node in returned_nodes],
                "edges": [{"from": g.vertex_properties["graphId"][source], "to": g.vertex_properties["graphId"][target]} for source, target in returned_edges]}
    node_types = {g.vertex_properties["graphId"][node]: g.vertex_properties["type"][node] for node in returned_nodes}
    is_seed = {g.vertex_properties["graphId"][node]: node in set(seed_graph_ids) for node in returned_nodes}
    returned_scores = {g.vertex_properties["graphId"][node]: None for node in returned_nodes}
    is_result = {g.vertex_properties["graphId"][node]: node in returned_nodes for node in returned_nodes}
    db_degrees = {g.vertex_properties["graphId"][node]: degrees[g.vertex_properties["graphId"][node]] for node in returned_nodes} 

    for node, score in best_drugs:
        returned_scores[g.vertex_properties["graphId"][node]] = score
    task_hook.set_results({
        "network": subgraph,
        "node_attributes": {
            "node_types": node_types, 
            "is_seed": is_seed, 
            "scores": returned_scores, 
            "is_result": is_result,
            "db_degrees": db_degrees,
            },
        # "traces_degree": traces_degree,
        "cancer_types": cancer_types,
        "cancer_dataset": cancer_dataset,
        "gene_interaction_dataset": gene_interaction_datasets[0],
        "drug_interaction_dataset": drug_interaction_datasets[0]
    })
