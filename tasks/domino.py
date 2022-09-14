from tkinter import E
from tasks.util.read_graph_tool_graph import read_graph_tool_graph
from tasks.task_hook import TaskHook
import sys
import itertools
import os.path
import os
import subprocess
import pathlib
import shutil
from pyensembl import EnsemblRelease


def symbols_to_ensembles(symbols):
    ensembl = EnsemblRelease(106)
    data = []
    mapping = {}
    for symbol in symbols:
        try:
            ensembles = ensembl.gene_ids_of_gene_name(symbol)
            data.extend(ensembles)
            for ensg in ensembles:
                mapping[ensg] = symbol
        except Exception as e:
            pass
    return data, mapping

def ensembles_to_symbols(ensembles, mapping):
    ensembl = EnsemblRelease(106)
    data = []
    for gene in ensembles:
        if gene in mapping:
            data.append(mapping[gene])
        else:
            try:
                 data.append(ensembl.gene_by_id(gene).to_dict()['gene_name'])
            except:
                pass
    return data

def write_gene_file(symbols, path):
    ensembles, ensg_to_symbol_mapping = symbols_to_ensembles(symbols)
    with open(path, 'a') as f:
        for ensg in ensembles:
            f.write(f'{ensg}\n')
    return ensg_to_symbol_mapping

def write_network_file(g, path):
    edges = []
    ensembl = EnsemblRelease(106)
    for edge in g.edges():
        name_a = g.vertex_properties["name"][edge.source()]
        name_b = g.vertex_properties["name"][edge.target()]
        
        try:
            ensemlbles_a = ensembl.gene_ids_of_gene_name(name_a)
            ensemlbles_b = ensembl.gene_ids_of_gene_name(name_b)
            
            for ensemble_a in ensemlbles_a:
                for ensemble_b in ensemlbles_b:
                    edges.append([ensemble_a, ensemble_b])
        except Exception as e:
            pass

    with open(path, 'a') as f:
        f.write('ID_interactor_A\tppi\tID_interactor_B\n')
        for a, b in edges:
            f.write(f'{a}\tppi\t{b}\n')

def read_output(output_path, symbol_id_map, g, ensg_to_symbol_mapping):

    f = open(output_path, "r")
    res = f.read()
    # string to lists
    modules = [x.strip('][').split(', ') for x in res.split('\n')]
    # ENSG to symbol
    modules = [ensembles_to_symbols(x, ensg_to_symbol_mapping) for x in modules]
    # to IDs
    modules = [[symbol_id_map[x] for x in l if x in symbol_id_map] for l in modules]

    # find edges in the modules
    edges = []
    for module in modules:
        for node1, node2 in itertools.product(module, repeat=2):
            if g.edge(node1, node2, add_missing=False) is not None:
                edges.append({'from': g.vertex_properties["graphId"][node1], 'to': g.vertex_properties["graphId"][node2]})
    return modules, edges

def domino_task(task_hook: TaskHook, token: str):

    seeds = task_hook.parameters["seeds"]
    seeds.sort()

    target = 'drug-target'

    cancer_types = task_hook.parameters.get("cancer_types", [])

    cancer_dataset = task_hook.parameters.get("cancer_dataset", "NCG6")

    gene_interaction_datasets = task_hook.parameters.get("gene_interaction_datasets", ["BioGRID"])
    
    max_deg = task_hook.parameters.get("max_deg", sys.maxsize)

    task_hook.set_progress(0 / 5.0, "Preparing input.")
    file_path = os.path.join(
        task_hook.data_directory,
        f"internal_{cancer_dataset}.gt"
    )

    g, seed_ids, cancer_node_ids, _, degrees = read_graph_tool_graph(
        file_path=file_path,
        gene_datasets=gene_interaction_datasets,
        drug_datasets=[],
        seeds=seeds,
        cancer_types=cancer_types,
        ignored_edge_types=[],
        max_deg=max_deg,
        target=target
    )

    symbols = [g.vertex_properties["name"][node] for node in seed_ids]

    task_hook.set_progress(1 / 5.0, "Writing genes.")

    folder = f'/tmp/{token}'
    pathlib.Path(folder).mkdir(parents=True, exist_ok=True) 
    
    gene_file_path = f'{folder}/{token}.txt'
    ensg_to_symbol_mapping = write_gene_file(symbols, gene_file_path)

    pathlib.Path( task_hook.data_directory+'domino').mkdir(parents=True, exist_ok=True) 
    network_path = os.path.join(
        task_hook.data_directory,
        'domino',
        f'internal_domino_{gene_interaction_datasets[0]}.sif'
        )
    
    slices_path = os.path.join(task_hook.data_directory, 'domino', f'internal_domino_{gene_interaction_datasets[0]}_sliced.txt')

    if not os.path.isfile(network_path):
        write_network_file(g, network_path)

    task_hook.set_progress(2 / 5.0, "Partitioning network.")
    if not os.path.isfile(slices_path):
        cmd = [
            'slicer', 
            '--network_file', network_path, 
            '--output_file', slices_path]
        subprocess.Popen(cmd).wait()
    
    output_folder = f'{folder}/{token}_output_folder'
    task_hook.set_progress(3 / 5.0, "Executing DOMINO.")
    cmd = [
        'domino', 
        '--active_genes_files', gene_file_path, 
        '--network_file', network_path, 
        '--slices_file', slices_path, 
        '--output_folder', output_folder]
    subprocess.Popen(cmd).wait()

    task_hook.set_progress(4 / 5.0, "Saving result.")

    out_file = os.path.join(output_folder, token, 'modules.out')
    symbol_id_map = {g.vertex_properties["name"][node]: node for node in seed_ids}
    modules, edges = read_output(out_file, symbol_id_map, g, ensg_to_symbol_mapping)

    node_ids = [node for module in modules for node in module]
    all_nodes = list(set(seed_ids + node_ids))
    subgraph = {"nodes": [g.vertex_properties["graphId"][node] for node in all_nodes],
                "edges": edges}
    node_types = {g.vertex_properties["graphId"][node]: g.vertex_properties["type"][node] for node in all_nodes}
    is_seed = {g.vertex_properties["graphId"][node]: node in node_ids for node in all_nodes}
    is_result = {g.vertex_properties["graphId"][node]: node in node_ids for node in all_nodes}
    db_degrees = {g.vertex_properties["graphId"][node]: degrees[g.vertex_properties["graphId"][node]] for node in all_nodes}
    cluster = {g.vertex_properties["graphId"][node]: i+1 for i, module in enumerate(modules) for node in module}

    task_hook.set_results({
        "network": subgraph,
        "node_attributes": {"node_types": node_types, "is_seed": is_seed, "is_result": is_result, "db_degrees": db_degrees, "cluster": cluster},
        "cancer_types": cancer_types,
        'cancer_dataset': cancer_dataset,
        'gene_interaction_dataset': gene_interaction_datasets[0],
        'drug_interaction_dataset': 'BioGRID'
    })

    # remove domino files
    shutil.rmtree(folder)

    task_hook.set_progress(5 / 5.0, "Done.")



