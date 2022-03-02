import graph_tool as gt
from graph_tool import topology
import pandas as pd
import glob
import os

from django.core.management import BaseCommand
import django

django.setup()


class Command(BaseCommand):
    def add_arguments(self, parser):
        pass

    def handle(self, *args, **kwargs):
        graph_file_paths = glob.glob('data/Common/internal_*.gt')

        for file_path in graph_file_paths:
            """
            Needs the .gt files, if they are not available or outdated/invalid, run "make_graphs.py" before this script
            """

            self.make_shortest_distance_file(file_path)

    def make_shortest_distance_file(self, file_path):

        print(file_path)
        output_filename = file_path[:-3] + '_shortest_distances_to_cancer_gene.csv'

        if os.path.isfile(output_filename):
            # check if file already exists
            return

        # load graph
        g = gt.load_graph(file_path)

        # drop duplicates due to cancer and non-cancer nodes
        to_remove = []
        cancer_nodes_graph_ids = []
        for node_id in range(g.num_vertices()):
            if g.vertex_properties['type'][node_id] == 'CancerNode':
                cancer_nodes_graph_ids.append(g.vertex_properties['graphId'][node_id])

        for node_id in range(g.num_vertices()):
            if g.vertex_properties['graphId'][node_id] in cancer_nodes_graph_ids \
                    and g.vertex_properties['type'][node_id] == 'Node':
                to_remove.append(node_id)

        g.remove_vertex(to_remove, fast=True)

        # get column and index names
        names = []
        for node_id in range(g.num_vertices()):
            names.append(g.vertex_properties['graphId'][node_id])

        # get distance iterator
        distances = topology.shortest_distance(g)

        # create distance matrix
        m = []
        for vector in distances:
            m.append(list(vector))
        # free space
        distances = None

        cancer_types = {}
        for node_id in range(g.num_vertices()):
            name = g.vertex_properties['graphId'][node_id]
            cancer_type = g.vertex_properties['cancer'][node_id]
            cancer_types[name] = cancer_type

        # put data in pandas DataFrame for more functionality
        df = pd.DataFrame(m, index=names, columns=names)

        # free space
        m = None

        # filter out non-cancer nodes
        cancer_nodes_map = []
        for node_id in range(g.num_vertices()):
            cancer_nodes_map.append(g.vertex_properties['type'][node_id] == 'CancerNode')
        df = df.iloc[cancer_nodes_map]

        # max val is all edges * all vertices in graph when no path found, filter out
        max_val = df.max().max()
        # collect shortest distances
        shortest_distances = []
        for name in names:
            # get min value without 0
            min_val = df[name][df[name] > 0].min()
            df_shortest_dist = df[name][(df[name] == min_val) & (df[name] != max_val)]
            for index, row in df_shortest_dist.to_frame().iterrows():
                for cancer_type in cancer_types[index].split(','):
                    shortest_distances.append(
                        {'source': name, 'target': index, 'dist': row[name], 'cancer_type': cancer_type}
                    )

        # free space
        df = None
        df_shortest_dist = None

        df_shortest_distances = pd.DataFrame(shortest_distances)
        df_shortest_distances.to_csv(output_filename)
