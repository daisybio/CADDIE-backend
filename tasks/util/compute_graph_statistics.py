from tasks.util.read_graph_tool_graph import read_graph_tool_graph
import graph_tool.topology as gtt
import sys
import pandas as pd
import argparse

def compute_graph_statistics(input_file_path, datasets, ignored_edge_types, output_file_path):
    g, _, bait_ids, _ = read_graph_tool_graph(input_file_path, [], datasets, ignored_edge_types, sys.maxsize)
    maxdeg = max([g.vertex(node).out_degree() for node in range(g.num_vertices())])
    non_bait_ids = [node for node in range(g.num_vertices()) if not node in set(bait_ids)]
    distances = {bait : gtt.shortest_distance(g, bait) for bait in bait_ids}
    distances_to_closest_baits = {non_bait : 1000000 for non_bait in non_bait_ids}
    for non_bait in non_bait_ids:
        for bait in bait_ids:
            distances_to_closest_baits[non_bait] = min(distances_to_closest_baits[non_bait], distances[bait][non_bait])
    closest_bait_names = {non_bait : [] for non_bait in non_bait_ids}
    for non_bait in non_bait_ids:
        for bait in bait_ids:
            if distances_to_closest_baits[non_bait] == distances[bait][non_bait]:
                closest_bait_names[non_bait].append(g.vertex_properties["name"][bait])
    data = {
        "UNIPROT" : [g.vertex_properties["name"][non_bait] for non_bait in non_bait_ids],
        "dist_to_closest_baits" : [distances_to_closest_baits[non_bait] for non_bait in non_bait_ids],
        "closest_baits" : [";".join(closest_bait_names[non_bait]) for non_bait in non_bait_ids]    
        }
    df = pd.DataFrame(data=data)
    df.to_csv(output_file_path, index=False)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser("Tool to compute distances to baits and maximum degree.")
    parser.add_argument("--input-file-path", required=True, help="Path to .gt or .graphml file.")
    parser.add_argument("--output-file-path", required=True, help="Path to .csv file where the distances to teh baits should be saved.")
    parser.add_argument("--datasets", nargs="+", default=[], help="Datasets that should be considered. If empty, all datasets are used.")
    parser.add_argument("--ignored-edge_types", nargs="+", default=[], help="Edge types that shoudl be ignored. If empty, all edge types are used.")
    args = parser.parse_args()
    compute_graph_statistics(args.input_file_path, args.datasets, args.ignored_edge_types, args.output_file_path)