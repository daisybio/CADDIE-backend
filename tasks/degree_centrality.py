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
    
    Parameters
    ----------
    strain_or_drugs : str
      The virus strain for which the analysis should be run, or the string literal "drugs"
      (if used for ranking drugs).
      
    seeds : list of str
      A list of UNIPROT IDs identifying the seed proteins.
      
    datasets : list of str
      List of datasets whose nodes and returned_edges should be considered for the analysis.
      If empty, all available datasets are employed.
      
    ignored_edge_types : list of str, optional (default: [])
      Edges whose types are contained in this list are ignored in the analysis.
      
    include_non_approved_drugs : bool, optional (default: False)
      If True, also non-approved drugs are considered when ranking drugs.
      
    result_size : int, optional (default: 20)
      The number of new candidate proteins to be returned.
      
    num_threads : int, optional (default: 1)
      Number of threads. Requires that graph_tool is compiled with OpenMP support.
      Should not be exposed in the frontend.
      
    Returns
    -------
    results : {"networks": list of dict, "node_attributes": list of dict}
      "networks": A one-element list containing the subgraph induced by the result 
        proteins and the seeds, along with the maximal score of all nodes, the maximal
        score off all nodes except virus proteins, and the maximal scores of all nodes
        except seed nodes.
      "node_attributes": A one-element list containing a dictionary with the following 
        attributes for all nodes in the returned network:
        "node_types": The type of the nodes (either , "host", or "drug").
        "is_seed": A flag that specifies whether the node is a seed.
        "scores": The un-normalized scores for all non-seed nodes (nan for the virus proteins).
      
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
    
    # Type: list of str
    # Semantics: Names of the seed proteins. Use UNIPROT IDs for host proteins, and 
    #            names of the for SARS_CoV2_<IDENTIFIER> (e.g., SARS_CoV2_ORF6) for
    #            virus proteins.
    # Reasonable default: None, has to be selected by user via "select for analysis"
    #            utility in frontend.
    # Acceptable values: UNIPROT IDs, identifiers of viral proteins.
    seeds = task_hook.parameters["seeds"]

    # Type: str
    # Semantics: sets the target of algorithms

    # Reasonable nodes, usually we target the nodes

    # Acceptable values: "nodes" or "drugs"
    target = task_hook.parameters.get("target", "drug-target")

    # Type: str.
    # Semantics: The virus strain for which the analysis should be run.
    # Example: "SARS_CoV2"
    # Reasonable default: None, has to be specified by the caller.
    # Acceptable values: "1", ...
    cancer_types = task_hook.parameters.get("cancer_types", [])
    
    # Type: str.
    # Semantics: The dataset which should be considered for the analysis.
    # Example: "1".
    # Reasonable default: "1".
    # Acceptable values: "1", "2"
    cancer_dataset = task_hook.parameters.get("cancer_dataset", "NCG6")

    # Type: str.
    # Semantics: The dataset which should be considered for the analysis.
    # Example: "1".
    # Reasonable default: "1".
    # Acceptable values: "1", "2"
    gene_interaction_datasets = task_hook.parameters.get("gene_interaction_datasets", ["BioGRID"])
    # gene_interaction_datasets = [models.InteractionGeneGeneDataset.objects.get(name__iexact=x) for x in gene_interaction_datasets]
    # gene_interaction_datasets = [f'{x.name}|{x.version}' for x in gene_interaction_datasets]

    # Type: str.
    # Semantics: The dataset which should be considered for the analysis.
    # Example: "1".
    # Reasonable default: "1".
    # Acceptable values: "1", "2"
    drug_interaction_datasets = task_hook.parameters.get("drug_interaction_datasets", ["BioGRID"])
    # drug_interaction_datasets = [models.InteractionGeneDrugDataset.objects.get(name__iexact=x) for x in drug_interaction_datasets]
    # drug_interaction_datasets = [f'{x.name}|{x.version}' for x in drug_interaction_datasets]
    
    # Type: list of str.
    # Semantics: Virus-host edge types which should be ignored for the analysis.
    # Example: ["Overexpression"].
    # Note: If empty, all available edge types are used. When ranking drugs, the 
    #       default [] should be used.
    # Reasonable default: [].
    # Acceptable values: "AP-MS", "overexpression".
    ignored_edge_types = task_hook.parameters.get("ignored_edge_types", [])
    
    # Type: bool
    # Semantics: Sepcifies whether should be included in the analysis when ranking drugs.
    # Example: False.
    # Reasonable default: False.
    # Has no effect unless trust_rank.py is used for ranking drugs.
    include_non_approved_drugs = task_hook.parameters.get("include_non_approved_drugs", False)
    
    # Type bool.
    # Semantics: Ignore viral proteins which are not selected as seeds.
    # Example: False.
    # Reasonable default: False.
    # Has no effect when the algorithm is used for ranking drugs.
    ignore_non_seed_baits = task_hook.parameters.get("ignore_non_seed_baits", False)
    
    # Type: int.
    # Semantics: All nodes with degree > max_deg * g.num_vertices() are ignored.
    # Example: 39.
    # Reasonable default: sys.maxsize.
    # Acceptable values: Positive integers.
    max_deg = task_hook.parameters.get("max_deg", sys.maxsize)
    
    # Type: int.
    # Semantics: Number of returned proteins.
    # Example: 20.
    # Reasonable default: 20.
    # Acceptable values: integers n with n > 0.
    result_size = task_hook.parameters.get("result_size", 20)
    
    # Type: int.
    # Semantics: Number of threads used for running the analysis.
    # Example: 1.
    # Reasonable default: 1.
    # Note: We probably do not want to expose this parameter to the user.
    num_threads = task_hook.parameters.get("num_threads", 1)

    # Type: boolean.
    # Semantics: Include nutraceutical drugs as candidates for drug search
    include_nutraceutical_drugs = task_hook.parameters.get("include_nutraceutical_drugs", False)

    include_only_ctrpv2_drugs = task_hook.parameters.get("include_only_ctrpv2_drugs", False)

    # Type: boolean.
    # Semantics: Include nutraceutical drugs as candidates for drug search
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