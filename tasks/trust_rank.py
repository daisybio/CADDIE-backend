from tasks.util.read_graph_tool_graph import read_graph_tool_graph
from tasks.util.scores_to_results import scores_to_results
from tasks.util.edge_weights import edge_weights
from tasks.task_hook import TaskHook
import graph_tool as gt
import graph_tool.centrality as gtc
import os.path
import sys
from caddie import models


def trust_rank(task_hook: TaskHook):
    r"""Computes TrustRank.
    
    The TrustRank [1_] is a node centrality measure which scores nodes in a network 
    based on how well they are connected to a (trusted) set of seed nodes (in our case: 
    proteins selected for analysis). It is a variant of the PageRank, where the node 
    centraities are initialized by assigning uniform probabilities to all seeds and zero 
    probabilities to all non-seed nodes.
    
    Parameters
    ----------
    strain_or_drugs : str
      The virus strain for which the analysis should be run, or the string literal "drugs"
      (if used for ranking drugs).
      
    seeds : list of str
      A list of UNIPROT IDs identifying the seed proteins.
      
    dataset : str
      List of dataset whose nodes and returned_edges should be considered for the analysis.
      If empty, all available dataset are employed.
      
    ignored_edge_types : list of str, optional (default: [])
      Edges whose types are contained in this list are ignored in the analysis.
      
    ignore_drug_non_seed_edges : bool, optional (default: False)
      If True, edges from non-seed host proteins to drugs are ignored when ranking drugs.
      
    include_non_approved_drugs : bool, optional (default: False)
      If True, also non-approved drugs are considered when ranking drugs.
      
    damping_factor : float, optional (default: 0.85)
      The employed damping factor (between 0 and 1). The larger the damping factor, 
      the faster the trust is propagated through the network.
      
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
    Let :math:`S` be the number of seeds, :math:`d` be the selected damping factor, 
    and :math:`N^-(u)` and :math:`deg^+(u)` be, respectively, the set of in-neighbors 
    and the out-degree of a node :math:`u`. Then the TrustRank :math:`TR(s)` of a 
    seed :math:`s` is recursively defined as
    
    .. math::
        TR(s) = \frac{1-d}{S} + d \sum_{u\in N^-(s)}\frac{TR(u)}{deg^+(u)}\text{,}
    
    while the TrustRank :math:`TR(v)` of a non-seed node $v$ is defined as follows: 
    
    .. math::
        TR(v) = d \sum_{u\in N^-(v)}\frac{TR(u)}{deg^+(u)} 
     
    The algorithm iteratively evaluates these equations until convergence. For undirected 
    graphs, the set of in-neighbours :math:`N^-(u)` and the out-degree :math:`deg^+(u)` 
    are substituted by the set of neighbors :math:`N(u)` and the degree :math:`deg(u)`, 
    respectively.
    
    This implementation is based on graph-tool, a very efficient Python package for network
    analysis with C++ backend and multi-threading support. Installation instructions for graph-tool
    can be found at https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions.
    
    References
    ----------
    .. [1] Z. Gyöngyi, H. Garcia-Molina, and J. O. Pedersen, Combating Web Spam with TrustRank,
       VLDB, 2004, pp. 576-587, https://doi.org/10.1016/B978-012088469-8.50052-8.
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

    # Acceptable values: "drug-target" or "drug"
    target = task_hook.parameters.get("target", "drug-target")

    # Type: str.
    # Semantics: The virus strain for which the analysis should be run.
    # Example: "SARS_CoV2"
    # Reasonable default: None, has to be specified by the caller.
    # Acceptable values: "1", ...
    cancer_types = task_hook.parameters.get("cancer_types", [])  
    # Type: list of str.
    # Semantics: The datasets which should be considered for the analysis.
    # Example: ["Krogan", "TUM"].
    # Note: If empty, all available datasets are used. When ranking drugs, the 
    #       default [] should be used.
    # Reasonable default: [].
    # Acceptable values: "Krogan", "TUM".
    cancer_dataset = task_hook.parameters.get("cancer_dataset", 'NCG6')

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
    
    # Type bool.
    # Semantics: Ignore viral proteins which are not selected as seeds.
    # Example: False.
    # Reasonable default: False.
    # Has no effect when the algorithm is used for ranking drugs.
    ignore_non_seed_baits = task_hook.parameters.get("ignore_non_seed_baits", False)
    
    # Type: bool
    # Semantics: Sepcifies whether also drugs targeting interactors of the seeds should be considered.
    # Example: False.
    # Reasonable default: False.
    # Has no effect unless trust_rank.py is used for ranking drugs.
    include_indirect_drugs = task_hook.parameters.get("include_indirect_drugs", False)
    
    # Type: bool
    # Semantics: Sepcifies whether should be included in the analysis when ranking drugs.
    # Example: False.
    # Reasonable default: False.
    # Has no effect unless trust_rank.py is used for ranking drugs.
    include_non_approved_drugs = task_hook.parameters.get("include_non_approved_drugs", False)
    
    # Type: double.
    # Semantics: Damping factor used by TrustRank. The larger the damping factor,
    #            the faster the trust is propagated through the network.
    # Example: 0.85.
    # Reasonable default: 0.85.
    # Acceptable values: doubles x with 0 < x < 1.
    damping_factor = task_hook.parameters.get("damping_factor", 0.85)
    
    # Type: int.
    # Semantics: Number of returned proteins.
    # Example: 20.
    # Reasonable default: 20.
    # Acceptable values: integers n with n > 0.
    result_size = task_hook.parameters.get("result_size", 20)
    
    # Type: int.
    # Semantics: All nodes with degree > max_deg * g.num_vertices() are ignored.
    # Example: 39.
    # Reasonable default: sys.maxsize.
    # Acceptable values: Positive integers.
    max_deg = task_hook.parameters.get("max_deg", sys.maxsize)

    # Type: float.
    # Semantics: Penalty parameter for hubs. Set edge weight to 1 + (hub_penalty / 2) (e.source.degree + e.target.degree)
    # Example: 0.5.
    # Reasonable default: 0.
    # Acceptable values: Floats between 0 and 1.
    hub_penalty = task_hook.parameters.get("hub_penalty", 0.0)
    
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
        raise ValueError(f'Could not find tissue.')

    drug_target_action = task_hook.parameters.get("drug_target_action", None)
    
    available_drugs = task_hook.parameters.get("available_drugs", None)
    if available_drugs is not None:
        available_drugs = [drug.lower() for drug in available_drugs]  
        
    # Parsing input file.
    task_hook.set_progress(0 / 4.0, "Parsing input.")
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
        include_nutraceutical_drugs=include_nutraceutical_drugs,
        only_atc_l_drugs=only_atc_l_drugs,
        target=target,
        drug_action=drug_target_action,
        available_drugs=available_drugs,
        include_only_ctrpv2_drugs=include_only_ctrpv2_drugs
    )

    task_hook.set_progress(1 / 4.0, "Computing edge weights.")
    weights = edge_weights(
        g,
        hub_penalty,
        mutation_cancer_type,
        expression_cancer_type,
        tissue,
        inverse=True,
    )
    
    # Set number of threads if OpenMP support is enabled.
    if gt.openmp_enabled():
        gt.openmp_set_num_threads(num_threads)
    
    # Call graph-tool to compute TrustRank.
    task_hook.set_progress(2 / 4.0, "Computing TrustRank.")
    trust = g.new_vertex_property("double")
    trust.a[seed_graph_ids] = 1.0 / len(seed_graph_ids)

    scores = gtc.pagerank(g, damping=damping_factor, pers=trust, weight=weights)

    # Compute and return the results.
    task_hook.set_progress(3 / 4.0, "Formating results.")
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
