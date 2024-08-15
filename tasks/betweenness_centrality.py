from tasks.util.read_graph_tool_graph import read_graph_tool_graph
from tasks.util.scores_to_results import scores_to_results
from tasks.util.edge_weights import edge_weights, custom_edge_weights
from tasks.task_hook import TaskHook
import graph_tool as gt
import graph_tool.topology as gtt
import os.path
import sys
from caddie import models


def betweenness_centrality(task_hook: TaskHook):
    r"""Computes betweenness centrality w.r.t. seed nodes.      
    Notes
    -----
    This implementation is based on graph-tool, a very efficient Python package for network
    analysis with C++ backend and multi-threading support. Installation instructions for graph-tool
    can be found at https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions.
    
    References
    ----------
    .. [1] L. Freeman, A Set of Measures of Centrality Based on Betweenness, Sociometry 40(1), 
       1977, pp. 35–41, https://doi.org/10.2307/3033543.
       [2] T. Kacprowski, N.T. Doncheva, M. Albrecht, NetworkPrioritizer: A Versatile Tool for 
       Network-Based Prioritization of Candidate Disease Genes or Other Molecules, Bioinformatics 29(11),
       2013, pp. 1471-1473, https://doi.org/10.1093/bioinformatics/btt164.  
    """
    
    seeds = task_hook.parameters["seeds"]

    target = task_hook.parameters.get("target", "drug-target")

    cancer_types = task_hook.parameters.get("cancer_types", [])

    cancer_dataset = task_hook.parameters.get("cancer_dataset", "NCG6")

    gene_interaction_datasets = task_hook.parameters.get("gene_interaction_datasets", ["BioGRID"])
    
    drug_interaction_datasets = task_hook.parameters.get("drug_interaction_datasets", ["BioGRID"])
    
    ignored_edge_types = task_hook.parameters.get("ignored_edge_types", [])
    
    include_indirect_drugs = task_hook.parameters.get("include_indirect_drugs", False)

    include_non_approved_drugs = task_hook.parameters.get("include_non_approved_drugs", False)

    ignore_non_seed_baits = task_hook.parameters.get("ignore_non_seed_baits", False)
    
    result_size = task_hook.parameters.get("result_size", 20)

    max_deg = task_hook.parameters.get("max_deg", sys.maxsize)

    hub_penalty = task_hook.parameters.get("hub_penalty", 0.0)

    num_threads = task_hook.parameters.get("num_threads", 1)

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

    task_hook.set_progress(0 / 3.0, "Parsing input.")
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
        ignored_edge_types=ignored_edge_types,
        max_deg=max_deg,
        ignore_non_seed_baits=ignore_non_seed_baits,
        include_indirect_drugs=include_indirect_drugs,
        include_non_approved_drugs=include_non_approved_drugs,
        target=target
    )

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

    # Set number of threads if OpenMP support is enabled.
    if gt.openmp_enabled():
        gt.openmp_set_num_threads(num_threads)

    # Call graph-tool to compute betweenness centrality.
    task_hook.set_progress(1 / 3.0, "Computing betweenness centralities.")
    scores = g.new_vertex_property("float")
    all_pairs = [(source, target) for source in seed_graph_ids for target in seed_graph_ids if source < target]
    for source, target in all_pairs:
        local_scores = g.new_vertex_property("float")
        num_paths = 0.0
        for path in gtt.all_shortest_paths(g, source, target, weights=weights):
            local_scores.a[path[1:-1]] += 1
            num_paths += 1
        if num_paths > 0:
            local_scores.a /= num_paths
            # normalize
            #local_scores.a = 2*local_scores.a / ((len(seeds)-1)*(len(seeds)-2))
        scores.a += local_scores.a

    # Compute and return the results.
    task_hook.set_progress(2 / 3.0, "Formating results.")

    task_hook.set_results(
        scores_to_results(
            target,
            result_size,
            g,
            seed_graph_ids,
            cancer_node_ids,
            drug_ids,
            scores,
            cancer_types,
            cancer_dataset,
            gene_interaction_datasets[0],
            drug_interaction_datasets[0],
            filter_paths, 
            degrees
        )
    )
