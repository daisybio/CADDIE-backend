import graph_tool as gt


def read_graph_tool_graph(file_path, gene_datasets, drug_datasets, seeds, cancer_types, ignored_edge_types, max_deg, node_key='graphId',
                          ignore_non_seed_baits=False, include_indirect_drugs=False, include_non_approved_drugs=False,
                          target='drug', include_nutraceutical_drugs=True, only_atc_l_drugs=False, drug_action=None,
                          available_drugs=None):
    r"""Reads a graph-tool graph from file.

    Dataset is given in file_path, hence we dont need to set it explicitly via a parameter

    Reads a graph-tool graph from graphml or gt file and returns is along
    with the internal IDs of the seed and viral seeds and the drugs.

    Also does preprocessing of the graph, e.g. filtering out unwanted nodes of unselected cancer types

    Parameters
    ----------
    file_path : str
      A string specifying the path to a graphml or gt file.

    seeds : list of str
      A list of entrez ids or uniprot acs identifying seeds
      TODO: backendIds cannot be used here in order to be open to user-input

    cancer_types : list of str
      A list of cancer type ids

    ignored_edge_types : list of str
      Edges whose types are contained in this list are ignored in the analysis., like "node-node"

    ignore_non_seed_baits : bool
      If True, cancer driver genes which are non selected as seeds are ignored.

    include_indirect_drugs : bool
      If True, edges from non-seed genes to drugs are ignored when ranking drugs.

    include_non_approved_drugs : bool
      If True, also non-approved drugs are included in the analysis

    Returns
    -------
    g : graph_tool.Graph
      The constructed graph.

    seed_ids : list of int
      The internal IDs of the seed seeds.

    cancer_driver_genes_ids : list of int
      The internal IDs of the cancer driver genes seeds.

    drug_ids : list of int
      The internal IDs of the drugs.
    """
    # Read the graph.
    g = gt.load_graph(file_path)
    
    cancer_types = list(map(str, cancer_types)) if cancer_types else []
    seeds_set = set(seeds)
    # Delete all nodes that are not contained in the selected cancer type
    deleted_nodes = []
    for node in range(g.num_vertices()):
        # we delete node if name is not in seeds and degree is bigger than max_degree
        # For undirected graphs, the “out-degree” is synonym for degree,
        # and in this case the in-degree of a vertex is always zero.
        if not g.vertex_properties[node_key][node] in seeds_set and g.vertex(node).out_degree() > max_deg:
            deleted_nodes.append(node)

        # node is CancerNode
        elif g.vertex_properties["type"][node] == "CancerNode":
            # delete nodes if ignore_non_seed_baits is True and cancerNode is not in seeds
            if ignore_non_seed_baits and g.vertex_properties[node_key][node] not in seeds_set:
                deleted_nodes.append(node)

            # remove cancerNodes if their cancer type is not selected
            elif not any(
                    [t in g.vertex_properties["cancer"][node].split(",") for t in cancer_types]):
                deleted_nodes.append(node)

            # for each cancer gene, there is also a normal gene. in case a gene is a cancer gene for any of the selected
            # cancer types, we see it as a cancer gene and remove equivalent in the normal genes.
            # otherwise, the cancer gene will be removed and the normal gene will be kept (see if statement above)
            else:
                # TODO: save in list and iterate over g a second time instead of always
                for node_j in range(g.num_vertices()):
                    if g.vertex_properties["type"][node_j] == "Node" and \
                            g.vertex_properties["graphId"][node] == g.vertex_properties["graphId"][node_j]:
                        deleted_nodes.append(node_j)
                        # we can break here since backend IDs in non-cancer genes are unique
                        break

        # remove all drugs from graph if we are not looking for drugs
        # in covex this is implemented via a separate file
        elif target != 'drug' and g.vertex_properties["type"][node] == "Drug":
            deleted_nodes.append(node)

        elif g.vertex_properties["type"][node] == "Drug":
          
            found = ['Tretinoin', 'Betulinic Acid', 'Thalidomide', 'Gossypol', 'Chlorambucil', 'Fluorouracil', 'Cimetidine', 'Azacitidine', 'Trifluoperazine', 'Paclitaxel', 'Tamoxifen', 'Teniposide', 'Sildenafil', 'Simvastatin', 'Procarbazine', 'Curcumin', 'Ciclopirox', 'Myricetin', 'Methotrexate', 'Lovastatin', 'Valdecoxib', 'Cyclophosphamide', 'Dacarbazine', 'Niclosamide', 'Prochlorperazine', 'Ifosfamide', 'Doxorubicin', 'Ouabain', 'Topotecan', 'Etoposide', 'Tanespimycin', '(S)-blebbistatin', 'Mitomycin', 'Tacrolimus', 'Staurosporine', 'Cerulenin', 'Dasatinib', 'Gefitinib', 'Erlotinib', 'Tacedinaline', 'Ciclosporin', 'Vorinostat', 'Sirolimus', 'Sitagliptin', 'Entinostat', 'Belinostat', 'Vincristine', 'Itraconazole', 'BIIB021', 'MMI-175', 'Dexamethasone', 'Imatinib', 'Decitabine', 'Axitinib', 'Selumetinib', 'Vandetanib', 'Sorafenib', 'Temozolomide', 'Navitoclax', 'Bexarotene', 'Nilotinib', 'Sunitinib', 'Bendamustine', 'Omacetaxine mepesuccinate', 'Oxaliplatin', 'Brefeldin A', 'Triptolide', 'Bortezomib', 'Zebularine', 'Pevonedistat', 'PAC-1', 'Gemcitabine', 'Olaparib', 'Indisulam', 'Tamatinib', 'BMS-754807', 'Tosedostat', 'Neratinib', 'Daporinad', 'MK-1775', 'Tivozanib', 'Tandutinib', 'Fulvestrant', 'Bafilomycin A1', 'Barasertib', 'Veliparib', 'Saracatinib', 'Afatinib', 'Cediranib', 'Canertinib', 'Obatoclax', 'Masitinib', 'Brivanib', 'SNS-032', 'PX-12', 'Nintedanib', 'Crizotinib', 'Foretinib', 'Regorafenib', 'OSI-930', 'Lenvatinib', '4-(6-HYDROXY-BENZO[D]ISOXAZOL-3-YL)BENZENE-1,3-DIOL', 'Lapatinib', 'SGX-523', 'Alisertib', 'Ruxolitinib', 'PHA-793887', 'Quizartinib', 'Fingolimod', 'Pazopanib', 'AZD-9684', 'Vorapaxar', 'Clofarabine', 'XL765', 'Bardoxolone methyl', 'Elocalcitol', 'Fluvastatin', 'Hyperforin', 'MK-0752', 'Semagacestat', 'Bosutinib', 'Temsirolimus', 'Abiraterone', 'Docetaxel', 'Nelarabine', 'Rigosertib', 'Silmitasertib', 'Momelotinib', 'Birinapant', 'Cabozantinib', 'Ibrutinib', 'Sotrastaurin', 'Tivantinib', 'Trametinib', 'Dinaciclib', 'Istradefylline', 'Alvocidib', 'Dabrafenib', 'Linsitinib']
            if g.vertex_properties["name"][node] not in found:
                deleted_nodes.append(node)
                continue

            if available_drugs is not None and g.vertex_properties["name"][node].lower() not in available_drugs:
                deleted_nodes.append(node)

            elif only_atc_l_drugs and not g.vertex_properties["is_antineoplastic_and_immunomodulating_agent"][node]:
                deleted_nodes.append(node)

            # remove nutraceutical drugs if include_nutraceutical is false
            elif include_nutraceutical_drugs and g.vertex_properties["is_nutraceutical"][node]:
                deleted_nodes.append(node)

    g.remove_vertex(deleted_nodes, fast=True)
    # Retrieve internal IDs of seed_ids and cancer_driver_gene_ids.
    ignored_edge_types = set(ignored_edge_types)
    seed_ids = []
    cancer_node_ids = []
    drug_ids = []
    # is matched is dict with nodes as keys and a boolean as values
    is_matched = {node: False for node in seeds_set}

    for node in range(g.num_vertices()):
        # if node is seed node, add to seed_ids and set matched to True
        if g.vertex_properties[node_key][node] in seeds_set:
            seed_ids.append(node)
            is_matched[g.vertex_properties[node_key][node]] = True

        node_type = g.vertex_properties["type"][node]
        # if node is cancer node, add to cancer_node_ids
        if node_type == "CancerNode":
            cancer_node_ids.append(node)
        # if node is Drug
        elif node_type == "Drug":
            status = g.vertex_properties["status"][node]
            if status == 'approved' or (include_non_approved_drugs and status == "unapproved"):
                drug_ids.append(node)

    # Check that all seeds have been matched and throw error, otherwise
    for node, found in is_matched.items():
        if not found:
            raise ValueError("Invalid seed node {}. No node named {} in {}.".format(node, node, file_path))


    # Delete edges that should be ignored or are not contained in the selected dataset.
    deleted_edges = []
    for edge in g.edges():
        if target != 'drug':
            # remove all edges to drugs
            if g.edge_properties["type"][edge] == "drug-gene":
              deleted_edges.append(edge)
              continue   

        if g.edge_properties["type"][edge] == "drug-gene":
          if not any([x in drug_datasets for x in g.edge_properties['dataset_name'][edge].split(',')]):
              deleted_edges.append(edge)
              continue
        else:
          if not any([x in gene_datasets for x in g.edge_properties['dataset_name'][edge].split(',')]):
                deleted_edges.append(edge)
                continue
        
        if drug_action != None and g.edge_properties["type"][edge] == "drug-gene":
            if drug_action.startswith('not'):
                if g.edge_properties['action'][edge] != drug_action[4:]:
                    continue
            else:
                if g.edge_properties['action'][edge] == drug_action:
                    continue
            deleted_edges.append(edge)
            continue

        # TODO is this wanted?
        # elif g.edge_properties["type"][edge] in ignored_edge_types:
        #     deleted_edges.append(edge)

        # check if cancer types are not related to edge
        # no cancer type in edge means that edge is not related to any cancer node in this dataset
        # TODO check if this is necessary, nodes are already removed
        if cancer_types and g.edge_properties["cancer"][edge] and \
                not any([t in g.edge_properties["cancer"][edge].split(",") for t in cancer_types]):
            # deleted_edges.append(edge)
            continue

    g.set_fast_edge_removal(fast=True)
    for edge in deleted_edges:
        g.remove_edge(edge)

    # save degree information before indirect edges to drugs are removed
    degrees = {g.vertex_properties["graphId"][node]: g.vertex(node).out_degree() for node in range(g.num_vertices())}

    if not include_indirect_drugs:
      deleted_edges = []
      for edge in g.edges():
        # if target is drug and edge type is drug-gene and indirect_drugs is false,
        # make sure that edges to indirect drugs are filtered out
        if g.edge_properties["type"][edge] == "drug-gene" and \
          (edge.target() not in seed_ids and edge.source() not in seed_ids):
            deleted_edges.append(edge)
      for edge in deleted_edges:
          g.remove_edge(edge)
      g.set_fast_edge_removal(fast=False)


    # Return the graph and the internal IDs of the seed_ids, cancer ids and drug ids. In addition, return the node degrees.
    return g, seed_ids, cancer_node_ids, drug_ids, degrees
