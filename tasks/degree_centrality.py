from tasks.util.read_graph_tool_graph import read_graph_tool_graph
from tasks.util.scores_to_results import scores_to_results
from tasks.task_hook import TaskHook
import graph_tool as gt
import os.path
import sys
from caddie import models


def degree_centrality(task_hook: TaskHook):
    r"""Computes degree centrality w.r.t. seed nodes.
    
    The degree centrality of a node :math:`u` in a graph :math:`G=(V,E)` is defined 
    as its degree :math:`deg(u)=|\{v\in V | (u,v)\in E\}|`. We here use the modified 
    version :math:`deg_S(u)=|\{v\in S | (u,v)\in E\}|` suggested in [2_], 
    where :math:`S\subseteq V` is a set of selected seed nodes.
    
    Notes
    -----
    This implementation is based on graph-tool, a very efficient Python package for network
    analysis with C++ backend and multi-threading support. Installation instructions for graph-tool
    can be found at https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions.
    
    References
    ----------
    .. [1] G. Sabidussi, The Centrality Index of a Graph, Psychometrika 31(4), 1966, pp. 581–603,
       https://doi.org/10.1007/bf02289527.
       [2] T. Kacprowski, N.T. Doncheva, M. Albrecht, NetworkPrioritizer: A Versatile Tool for 
       Network-Based Prioritization of Candidate Disease Genes or Other Molecules, Bioinformatics 29(11),
       pp. 1471-1473, https://doi.org/10.1093/bioinformatics/btt164.  
    """
    
    seeds = task_hook.parameters["seeds"]

    target = task_hook.parameters.get("target", "drug-target")

    cancer_types = task_hook.parameters.get("cancer_types", [])
    
    cancer_dataset = task_hook.parameters.get("cancer_dataset", "NCG6")

    gene_interaction_datasets = task_hook.parameters.get("gene_interaction_datasets", ["BioGRID"])
   
    drug_interaction_datasets = task_hook.parameters.get("drug_interaction_datasets", ["BioGRID"])
    
    ignored_edge_types = task_hook.parameters.get("ignored_edge_types", [])
    
    include_non_approved_drugs = task_hook.parameters.get("include_non_approved_drugs", False)
    
    ignore_non_seed_baits = task_hook.parameters.get("ignore_non_seed_baits", False)
    
    max_deg = task_hook.parameters.get("max_deg", sys.maxsize)
    
    result_size = task_hook.parameters.get("result_size", 20)
    
    num_threads = task_hook.parameters.get("num_threads", 1)

    include_nutraceutical_drugs = task_hook.parameters.get("include_nutraceutical_drugs", False)

    include_only_ctrpv2_drugs = task_hook.parameters.get("include_only_ctrpv2_drugs", False)

    only_atc_l_drugs = task_hook.parameters.get("only_atc_l_drugs", False)

    filter_paths = task_hook.parameters.get("filter_paths", True)

    drug_target_action = task_hook.parameters.get("drug_target_action", None)
    
    available_drugs = task_hook.parameters.get("available_drugs", None)
    if available_drugs is not None:
        available_drugs = [drug.lower() for drug in available_drugs]
    
    # Parsing input file.
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
        include_indirect_drugs=False,
        include_non_approved_drugs=include_non_approved_drugs,
        include_nutraceutical_drugs=include_nutraceutical_drugs,
        only_atc_l_drugs=only_atc_l_drugs,
        target=target,
        drug_action=drug_target_action,
        available_drugs=available_drugs,
        include_only_ctrpv2_drugs=include_only_ctrpv2_drugs
    )
    
    # Set number of threads if OpenMP support is enabled.
    if gt.openmp_enabled():
        gt.openmp_set_num_threads(num_threads)
    
    # Call graph-tool to compute TrustRank.
    task_hook.set_progress(1 / 3.0, "Computing degree centralities.")
    scores = g.new_vertex_property("float")
    for node in seed_graph_ids:
        for nb in g.get_all_neighbors(node):
            scores.a[nb] += 1
    
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