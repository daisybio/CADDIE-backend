from scripts.minimum_steiner_tree.minimum_steiner_tree import min_steiner_tree
from tasks.util.read_graph_tool_graph import read_graph_tool_graph
from caddie import models, serializers
from django import db
from django.core.management import BaseCommand
import django
import multiprocessing

django.setup()

KERNEL = 6


def calculate_min_spanning_tree(params):
    print(params)
    dataset_object = params['cancer_dataset']
    gene_interaction_dataset = params['gene_interaction_dataset']
    drug_interaction_dataset = params['drug_interaction_dataset']
    cancer_type = params['cancer_type_id']

    file_path = f'data/networks/internal_{dataset_object.name}.gt'
    cancer_type_object = models.CancerType.objects.get(id=cancer_type)

    cancer_gene_objects = models.CancerGeneEntity.objects.filter(
        cancer_dataset_id=dataset_object,
        cancer_type_id=cancer_type_object
    )
    cancer_genes = serializers.CancerGeneEntitySerializer(many=True).to_representation(cancer_gene_objects)
    seeds = [g['graphId'] for g in cancer_genes]

    g, seed_ids, cancer_node_ids, drug_ids, degrees = read_graph_tool_graph(
        file_path=file_path,
        gene_datasets=[gene_interaction_dataset],
        drug_datasets=[drug_interaction_dataset],
        seeds=seeds,
        cancer_types=[cancer_type],
        ignored_edge_types=[],
        max_deg=99999999999,
        target='drug-target',
    )

    min_span_tree, vfilt = min_steiner_tree(g, seed_ids)  # returns a tree of type graph_tool.Graph

    min_spanning_tree_graph_ids = [
        g.vertex_properties['graphId'][node] for node in range(g.num_vertices()) if vfilt[node]
    ]

    # if cancer type has only one gene or genes do not have connections, they might not appear in result
    # here we add the lost genes manually to the result
    for graph_id in seeds:
        if graph_id not in min_spanning_tree_graph_ids:
            min_spanning_tree_graph_ids.append(graph_id)

    dataset_object = models.CancerDataset.objects.filter(name=dataset_object.name).first()
    interaction_dataset_object = models.InteractionGeneGeneDataset.objects.filter(
        name=gene_interaction_dataset).first()

    for graph_id in min_spanning_tree_graph_ids:
        gene = models.Gene.objects.get(id=graph_id[1:])
        models.MinSpanningTree.objects.get_or_create(
            cancer_type=cancer_type_object,
            gene=gene,
            cancer_dataset=dataset_object,
            edge_dataset=interaction_dataset_object,
            is_cancer=graph_id in seeds,
        )


class Command(BaseCommand):
    def add_arguments(self, parser):
        pass

    def handle(self, *args, **kwargs):
        cancer_datasets = models.CancerDataset.objects.all()
        # cancer_dataset_names = [(e.id, e.name) for e in cancer_datasets]

        cancer_gene_interactions = models.InteractionGeneGeneDataset.objects.all()
        cancer_gene_interaction_names = [e.name for e in cancer_gene_interactions]

        # does not matter, placeholder
        drug_interaction_dataset = 'BioGRID'

        parameters = []
        for cancer_dataset in cancer_datasets:
            cancer_types = models.RelationDatasetType.objects\
                .filter(cancer_dataset_id=cancer_dataset.id).order_by('cancer_type_id__name')
            cancer_type_ids = [e.cancer_type_id.id for e in cancer_types]
            for cancer_type_id in cancer_type_ids:
                cancer_type_id = str(cancer_type_id)
                for gene_interaction_dataset in cancer_gene_interaction_names:
                    parameters.append({
                        'cancer_dataset': cancer_dataset,
                        'gene_interaction_dataset': gene_interaction_dataset,
                        'drug_interaction_dataset': drug_interaction_dataset,
                        'cancer_type_id': cancer_type_id
                    })

        db.connections.close_all()
        pool = multiprocessing.Pool(KERNEL)
        pool.map(calculate_min_spanning_tree, parameters)
