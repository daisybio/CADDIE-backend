import json
import random
import string
import math
import csv

from django.db.models import Q
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.views import APIView
from rest_framework.exceptions import APIException
from django.http import HttpResponse, HttpResponseBadRequest

import networkx as nx

import caddie.models as models
import caddie.serializers as serializers

from caddie.tasks import start_task, refresh_from_redis, task_stats, task_result

from caddie.utils.vcfUtils import VCFUtils


def pointsInCircum(r, n=100):
    pi = math.pi
    return [(math.cos(2 * pi / n * x) * r, math.sin(2 * pi / n * x) * r) for x in range(0, n + 1)]


class CancerDatasetView(APIView):
    """
    List all the cancer datasets in the database
    """

    def get(self, request):
        """
        Return a list of all cancer datasets with their attributes
        """

        cancer_datasets_objects = models.CancerDataset.objects.all()
        cancer_datasets = serializers.CancerDatasetSerializer(many=True).to_representation(cancer_datasets_objects)
        # TODO store this information in summary table to speed up loading time
        # get count of cancer genes in this dataset
        for cancer_dataset in cancer_datasets:
            dataset_backendId = cancer_dataset['backendId']
            n = models.CancerGeneEntity.objects.filter(
                cancer_dataset_id=dataset_backendId
            ).distinct('cancer_gene_id').count()
            cancer_dataset['count'] = n

        return Response(cancer_datasets)


class InternalCancerGeneView(APIView):
    def get(self, request):
        cancer_dataset_id = json.loads(request.query_params.get('dataset'))
        cancer_gene_objects = models.CancerGeneEntity.objects.filter(
            cancer_dataset_id=int(cancer_dataset_id)
        )
        cancer_genes = serializers.CancerGeneEntitySerializer(many=True).to_representation(cancer_gene_objects)
        return Response(cancer_genes)


class InternalGeneInteractionView(APIView):
    def get(self, request):
        interaction_dataset_id = json.loads(request.query_params.get('dataset'))
        interaction_objects = models.GeneGeneInteraction.objects.filter(
            interaction_dataset_id=interaction_dataset_id
        )
        interactions = serializers.GeneGeneInteractionSerializer(many=True)\
            .to_representation(interaction_objects)
        return Response(interactions)


class InternalGeneDrugInteractionView(APIView):
    def get(self, request):
        interaction_dataset_id = json.loads(request.query_params.get('dataset'))
        interaction_objects = models.GeneDrugInteraction.objects.filter(
            interaction_dataset_id=interaction_dataset_id
        )
        interactions = serializers.GeneDrugInteractionSerializer(many=True)\
            .to_representation(interaction_objects)
        return Response(interactions)


class InteractionGeneDatasetView(APIView):
    """
    List all the gene interaction datasets in the database
    """

    def get(self, request):
        """
        Return a list of all interaction datasets with their attributes
        """

        interaction_datasets_objects = models.InteractionGeneGeneDataset.objects.all()
        interaction_datasets = serializers.InteractionDatasetSerializer(many=True)\
            .to_representation(interaction_datasets_objects)

        # TODO store this information in summary table to speed up loading time
        # get count of interactions in this dataset
        for interaction_dataset in interaction_datasets:
            dataset_backendId = interaction_dataset['backendId']
            n = models.GeneGeneInteraction.objects.filter(interaction_dataset_id=dataset_backendId).count()
            interaction_dataset['count'] = n

        return Response(interaction_datasets)


class InteractionDrugDatasetView(APIView):
    """
    List all the drug interaction datasets in the database
    """

    def get(self, request):
        """
        Return a list of all interaction datasets with their attributes
        """

        interaction_datasets_objects = models.InteractionGeneDrugDataset.objects.all()
        interaction_datasets = serializers.InteractionDatasetSerializer(many=True)\
            .to_representation(interaction_datasets_objects)

        # TODO store this information in summary table to speed up loading time
        # get count of interactions in this dataset
        for interaction_dataset in interaction_datasets:
            dataset_backendId = interaction_dataset['backendId']
            n = models.GeneDrugInteraction.objects.filter(interaction_dataset_id=dataset_backendId).count()
            interaction_dataset['count'] = n

        return Response(interaction_datasets)


class CancerTypeView(APIView):
    """
    Lists all cancer types
    """
    def get(self, request):
        """
        Return a list of all cancer types with their attributes
        """
        cancer_dataset_id = None
        if request.query_params.get('dataset'):
            cancer_dataset_id = json.loads(request.query_params.get('dataset'))

        cancer_types_objects = models.RelationDatasetType.objects.filter(cancer_dataset_id=cancer_dataset_id)\
            .order_by('cancer_type_id__name')
        cancer_types = serializers.RelationDatasetTypeSerializer(many=True).to_representation(cancer_types_objects)

        return Response(cancer_types)


class NodeInteractionLookup(APIView):
    """
    Returns all interactions for given backendId
    """

    def get(self, request):
        # TODO validate request in separate function
        if request.query_params.get('interactionDataset'):
            interaction_dataset = json.loads(request.query_params.get('interactionDataset'))
        else:
            return APIException('interactionDataset is missing')

        if request.query_params.get('backendId'):
            backend_id = json.loads(request.query_params.get('backendId'))
        else:
            return APIException('backendId is missing')

        interaction_objects = models.GeneGeneInteraction.objects.filter(
            (Q(gene_a=backend_id) | Q(gene_b=backend_id)) & Q(interaction_dataset_id=interaction_dataset)
        )

        interactions = serializers.GeneGeneInteractionSerializer(many=True).to_representation(interaction_objects)

        # TODO also return related genes?

        return Response({
            'interactions': interactions
        })


class RelatedCancerTypeView(APIView):

    def get(self, request):
        """
        Lookup for names of related cancer types for node in frontend (even if node is non-cancer driver gene)
        """
        if request.query_params.get('dataset'):
            cancer_dataset_id = json.loads(request.query_params.get('dataset'))
        else:
            return APIException('dataset is missing')

        if request.query_params.get('node_backendId'):
            node_backendId = json.loads(request.query_params.get('node_backendId'))
        else:
            return APIException('node_backendId is missing')

        cancer_gene_objects = models.CancerGeneEntity.objects\
            .filter(cancer_dataset_id=cancer_dataset_id, cancer_gene_id=node_backendId)

        cancer_genes = serializers.CancerGeneEntitySerializer(many=True).to_representation(cancer_gene_objects)

        # TODO create Serializer for that to speed up
        cancer_types = [{'name': gene['type'], 'backendId': gene['type_id']} for gene in cancer_genes]

        return Response({
            'cancer_types': cancer_types
        })


class ComorbidityNodeView(APIView):
    def get(self, request):
        """

        :param request:
        :return:
        """
        if request.query_params.get('node_backendId'):
            node_backendId = json.loads(request.query_params.get('node_backendId'))
        else:
            return APIException('node_backendId is missing')

        disease_gene_interaction_objects = models.DiseaseGeneInteractions.objects.filter(
            gene=node_backendId
        )

        disease_gene_interactions = serializers.DiseaseGeneInteractionsSerializer(many=True)\
            .to_representation(disease_gene_interaction_objects)
        return Response({
            'interactions': disease_gene_interactions
        })


class ComorbidityCancerTypeView(APIView):
    """
    Returns Data for horizontal barplot.
    Expects the keys "cancerDatasetBackendId" and "cancerTypeBackendIds" (list)
    """
    def get(self, request):

        if request.query_params.get('cancerDatasetBackendId'):
            cancerDatasetBackendId = json.loads(request.query_params.get('cancerDatasetBackendId'))
        else:
            return APIException('cancerDatasetBackendId is missing')

        if request.query_params.get('cancerTypeBackendIds'):
            cancerTypeBackendIds = json.loads(request.query_params.get('cancerTypeBackendIds')).split(',')
        else:
            return APIException('cancerTypeBackendIds is missing')

        # get genes of cancer type
        cancer_gene_objects = models.CancerGeneEntity.objects.filter(
            cancer_dataset_id=cancerDatasetBackendId,
            cancer_type_id__in=cancerTypeBackendIds
        )
        cancer_genes = serializers.CancerGeneEntitySerializer(many=True).to_representation(cancer_gene_objects)
        cancer_gene_ids = [gene['backendId'] for gene in cancer_genes]

        # get related diseases
        disease_gene_interaction_objects = models.DiseaseGeneInteractions.objects.filter(gene__in=cancer_gene_ids)
        disease_gene_interactions = serializers.DiseaseGeneInteractionsSerializer(many=True)\
            .to_representation(disease_gene_interaction_objects)

        # count occurances
        diseases_counts = {}
        disease_objects = {}
        for edge in disease_gene_interactions:
            if edge['disease_name'] in diseases_counts:
                diseases_counts[edge['disease_name']] += 1
            else:
                diseases_counts[edge['disease_name']] = 1
                disease_objects[edge['disease_name']] = serializers.DiseaseSerializer()\
                    .to_representation(models.Disease.objects.get(id=edge['disease_backendId']))

        return Response({
            'counts': diseases_counts,
            'n_cancer_genes': len(cancer_genes),
            'disease_objects': disease_objects
        })


class Network(APIView):
    """
    Collect and return data for protein protein interaction network of cancer driver genes
    """

    def get(self, request):
        """
        Fetches cancer genes and interactions from the database and returns all data relevant for the network
        :param request:
        :return: network data
        """

        # TODO validate request in separate function
        if request.query_params.get('dataset'):
            cancer_dataset_id = json.loads(request.query_params.get('dataset'))
        else:
            return APIException('dataset is missing')

        if request.query_params.get('interactionDataset'):
            interaction_dataset = json.loads(request.query_params.get('interactionDataset'))
        else:
            return APIException('interactionDataset is missing')

        if request.query_params.get('cancerTypes'):
            cancer_type_ids = json.loads(request.query_params.get('cancerTypes'))
            # cancer_type_ids is a stringify-ed list of cancer type ids
            cancer_type_ids_list = cancer_type_ids.split(',')
        else:
            return APIException('cancerTypes is missing')

        # look up nodes in pre-calculated min-spanning-tree table which are minimum spanning trees based on all
        # cancer genes in a cancer type, but every gene just once if they occur in multiple types
        node_objects = models.MinSpanningTree.objects.filter(
            cancer_dataset=cancer_dataset_id,
            edge_dataset=interaction_dataset,
            cancer_type__in=cancer_type_ids_list
        ).distinct('gene')
        # get all nodes from the node ids
        cancer_nodes = []
        nodes = []
        node_ids = []
        coords = pointsInCircum(2000, len(node_objects))

        for index, node in enumerate(node_objects):
            if node.is_cancer:
                gene_object = models.CancerGeneEntity.objects.filter(
                    cancer_gene_id=node.gene,
                    cancer_type_id__in=cancer_type_ids_list,
                    cancer_dataset_id=cancer_dataset_id
                ).first()
                gene = serializers.CancerGeneEntitySerializer().to_representation(gene_object)

                # TODO save this in db
                # for fixed positioning in graph
                gene['x'] = coords[index][0]
                gene['y'] = coords[index][1]
                cancer_nodes.append(gene)
            else:
                gene_object = models.Gene.objects.get(id=node.gene.id)
                gene = serializers.GeneSerializer().to_representation(gene_object)

                gene['x'] = coords[index][0]
                gene['y'] = coords[index][1]

                nodes.append(gene)
            node_ids.append(node.gene)

        # collect interactions between nodes
        gene_interaction_objects = models.GeneGeneInteraction.objects.filter(
            gene_a__in=node_ids,
            gene_b__in=node_ids,
            interaction_dataset_id=interaction_dataset
        )
        gene_interactions = serializers.GeneGeneInteractionSerializer(many=True)\
            .to_representation(gene_interaction_objects)

        nodes_total_n = models.Gene.objects.all().count() - len(cancer_nodes)

        network = {
            'cancerNodesSup': [],
            'cancerNodes': cancer_nodes,
            'nodes': nodes,
            'nodesSup': [],  # not used yet
            'nodesTotalN': nodes_total_n,
            'interactions': gene_interactions,
            'interactionsSup': [],  # not used yet
            'interactionsTotalN': len(gene_interactions)
        }

        return Response(network)


class VCFSeedLookup(APIView):

    def post(self, request):
        """
        Recieve a VCF file and
        :param request:
        :return:
        """
        if request.data.get('file_data'):
            file_data = json.loads(request.data.get('file_data'))
        else:
            return APIException('fileData is missing')

        if request.data.get('threshold'):
            threshold = json.loads(request.data.get('threshold'))
        else:
            return APIException('threshold is missing')

        # white spaces in possible annotation break the reader
        seeds_hugo = VCFUtils.filter_file(file_data, threshold)
        seed_objects = models.Gene.objects.filter(name__in=seeds_hugo)
        seeds = serializers.GeneSerializer(many=True).to_representation(seed_objects)
        return Response({'data': seeds})


class SummaryView(APIView):

    def get(self, request):
        ret = {
            'nGenes': models.Gene.objects.count(),
            'nCancerDriverGenes': models.CancerGeneEntity.objects.distinct('cancer_gene_id__entrez_id').count(),
            'nDrugs': models.Drug.objects.count(),
        }
        return Response(ret)


@api_view(['POST'])
def task_summarize(request):
    chars = string.ascii_lowercase + string.ascii_uppercase + string.digits
    token_str = ''.join(random.choice(chars) for _ in range(32))

    # parameters = request.data['parameters']

    # for source_task_token in parameters['source_tasks']:
    #     task = models.Task.objects.get(token=source_task_token)
    #     parameters = always_merger.merge(parameters, json.loads(task.parameters))

    task = models.Task.objects.create(
        token=token_str,
        target='summarize',
        algorithm='summary',
        parameters=json.dumps(request.data['parameters'])
    )
    start_task(task)
    task.save()

    return Response({
        'token': token_str,
    })


class TaskView(APIView):

    def post(self, request):
        # we take the cancer dataset name instead of id because we need to map it to the .gt files,
        # which are labeled by name in order to find the correct file even after reconstructing the db
        # which might shift indices
        chars = string.ascii_lowercase + string.ascii_uppercase + string.digits
        token_str = ''.join(random.choice(chars) for _ in range(32))

        parameters = request.data['parameters']

        cancer_type_names = parameters.get('cancer_type_names', [])
        cancer_type_objects = models.CancerType.objects.filter(name__in=cancer_type_names)
        parameters['cancer_types'] = [cancer_type_object.id for cancer_type_object in cancer_type_objects]

        # add target key to task parameters for simplicity
        parameters['target'] = request.data['target']
        task = models.Task.objects.create(
            token=token_str,
            target=request.data['target'],
            algorithm=request.data['algorithm'],
            parameters=json.dumps(parameters)
        )
        start_task(task)
        task.save()

        return Response({
            'token': token_str,
        })

    def get(self, request):
        token_str = request.query_params['token']
        task = models.Task.objects.get(token=token_str)

        if not task.done and not task.failed:
            refresh_from_redis(task)
            task.save()

        return Response({
            'token': task.token,
            'info': serializers.TaskSerializer().to_representation(task),
            'stats': task_stats(task),
        })


@api_view(['POST'])
def tasks_view(request):
    tokens = json.loads(request.data.get('tokens', '[]'))
    tasks = models.Task.objects.filter(token__in=tokens).order_by('-created_at').all()
    tasks_info = []
    for task in tasks:
        if not task.done and not task.failed:
            refresh_from_redis(task)
            task.save()

        tasks_info.append({
            'token': task.token,
            'info': serializers.TaskStatusSerializer().to_representation(task),
            'stats': task_stats(task),
        })
    return Response(tasks_info)


def infer_node_type_and_details(node, cancer_types):
    node_id = int(node[1:])
    if node[0] == 'g':
        # TODO differentiate between node and cancer node
        try:
            # take first element since filter returns list
            cancer_gene = models.CancerGeneEntity.objects.filter(
                cancer_gene_id__id=node_id,
                cancer_type_id__id__in=cancer_types
            ).first()
            # TODO: Add dataset to task and choose correct effect
            return 'CancerNode', cancer_gene
        except models.CancerGeneEntity.DoesNotExist:
            pass

        try:
            gene = models.Gene.objects.get(id=node_id)
            return 'Node', gene
        except models.Gene.DoesNotExist:
            pass

    elif node[0] == 'd':
        try:
            drug = models.DrugEntity.objects.get(drug_id=node_id)
            return 'Drug', drug
        except models.DrugEntity.DoesNotExist:
            pass

    return None, None


def get_closest_cancer_genes_for_gene(
        cancer_dataset,
        gene_gene_dataset,
        cancer_types,
        gene_id):
    """
    Finds closest cancer gene for given gene

    Distance is measured by nodes in between

    :param datasets:
    :param protein_ac:
    :return:
    """
    cancer_gene_objects = models.ShortestDistanceGeneToCancerGene.objects.filter(
        cancer_dataset__name=cancer_dataset,
        gene_interaction_dataset__name=gene_gene_dataset,
        cancer_type__in=cancer_types,
        gene_a=gene_id
    ).all()

    cancer_genes = serializers.ShortestDistanceCancerGeneSerializer(many=True).to_representation(cancer_gene_objects)
    return cancer_genes


def get_closest_cancer_genes_for_drug(
        cancer_dataset,
        drug_interaction_dataset,
        cancer_types,
        drug_id):
    """
    Finds closest cancer gene for given gene

    Distance is measured by nodes in between

    :param datasets:
    :param protein_ac:
    :return:
    """
    cancer_gene_objects = models.ShortestDistanceDrugToCancerGene.objects.filter(
        cancer_dataset__name=cancer_dataset,
        gene_drug_interaction_dataset__name=drug_interaction_dataset,
        cancer_type__in=cancer_types,
        drug=drug_id
    ).all()
    cancer_genes = serializers.ShortestDistanceDrugToCancerGeneSerializer(many=True)\
        .to_representation(cancer_gene_objects)
    return cancer_genes


def get_unique_bait_distances(baits):
    """
    Gets all unique baits in bait list

    :param baits:
    :return:
    """
    unique_bait_ids = set()
    unique_baits = []
    for b in baits:
        if b['gene_name'] not in unique_bait_ids:
            unique_bait_ids.add(b['gene_name'])
            unique_baits.append(b)
    return unique_baits


@api_view()
def result_view(request):
    view = request.query_params.get('view')
    fmt = request.query_params.get('fmt')
    token_str = request.query_params['token']
    task = models.Task.objects.get(token=token_str)
    result = task_result(task)

    cancer_types = result.get('cancer_types')
    cancer_types = list(map(int, cancer_types))

    cancer_dataset = result.get('cancer_dataset')
    gene_interaction_dataset = result.get('gene_interaction_dataset')
    drug_interaction_dataset = result.get('drug_interaction_dataset')

    node_attributes = result.get('node_attributes')

    if not node_attributes:
        node_attributes = {}
        result['node_attributes'] = node_attributes

    cancer_genes = []
    genes = []
    drugs = []

    network = result['network']

    node_types = node_attributes.get('node_types')
    if not node_types:
        node_types = {}
        node_attributes['node_types'] = node_types

    is_seed = node_attributes.get('is_seed')
    cluster = node_attributes.get('cluster', {})
    db_degrees = node_attributes.get('db_degrees')

    if not is_seed:
        is_seed = {}
        node_attributes['is_seed'] = is_seed

    scores = node_attributes.get('scores', {})

    node_details = {}
    node_attributes['details'] = node_details

    # parameters = json.loads(task.parameters)
    seeds = [k for k, v in node_attributes['is_seed'].items() if v]
    # nodes is list of gene names
    nodes = network['nodes']
    edges = network['edges']

    for node_graphId in nodes:
        node_backendId = int(node_graphId[1:])
        is_seed[node_graphId] = node_graphId in seeds
        node_type = node_types.get(node_graphId)
        pvd_entity = None
        details_s = None

        if not node_type:
            node_type, pvd_entity = infer_node_type_and_details(node_graphId, cancer_types)
        else:
            if node_type == 'Node':
                pvd_entities = [models.Gene.objects.get(id=node_backendId)]
            elif node_type == 'CancerNode':
                # can be more than one result due to cancer types
                pvd_entities = models.CancerGeneEntity.objects.filter(
                    cancer_dataset_id__name=cancer_dataset,
                    cancer_gene_id=node_backendId,
                    cancer_type_id__in=cancer_types
                )

            elif node_type == 'Drug':
                pvd_entities = [models.Drug.objects.get(id=node_backendId)]

        # we create 'pvd_entities' which is mostly just one element in a list for the case
        # that one cancer node occurs in multiple cancer types. In this case, we iterate over the x elements
        # for drugs and genes (and most cancer genes), the list will just contain one element
        for pvd_entity in pvd_entities:
            # this step is to check if no entry was found it db
            if not node_type or not pvd_entity:
                print('not found in db (should not occur) ' + node_graphId)
                continue

            if node_type == 'Node':
                details_s = serializers.GeneSerializer().to_representation(pvd_entity)
            elif node_type == 'CancerNode':
                details_s = serializers.CancerGeneEntitySerializer().to_representation(pvd_entity)
            elif node_type == 'Drug':
                details_s = serializers.DrugSerializer().to_representation(pvd_entity)

            if node_graphId not in node_types:
                node_types[node_graphId] = node_type
            else:
                # TODO add node type
                node_types[node_graphId] = node_type

            if scores.get(node_graphId) is not None:
                details_s['score'] = scores.get(node_graphId, None)

            details_s['cluster'] = cluster[node_graphId] if node_graphId in cluster else 0

            if node_type == 'Node':
                closest_baits = get_closest_cancer_genes_for_gene(
                    cancer_dataset,
                    gene_interaction_dataset,
                    cancer_types,
                    node_backendId
                )

                minimum_distance = None
                if len(closest_baits) > 0:
                    minimum_distance = min([cb['distance'] for cb in closest_baits])
                details_s['closest_cancer_genes'] = ','.join([cb['gene_b_name'] for cb in closest_baits])
                details_s['closest_distance'] = minimum_distance if minimum_distance else ''

            elif node_type == 'Drug':
                ds_closest_baits = []
                for edge in edges:
                    if edge['from'] == node_graphId:
                        ds_closest_baits.extend(get_closest_cancer_genes_for_drug(
                            cancer_dataset,
                            drug_interaction_dataset,
                            cancer_types,
                            int(edge['to'][1:])
                        ))
                    elif edge['to'] == node_graphId:
                        ds_closest_baits.extend(get_closest_cancer_genes_for_drug(
                            cancer_dataset,
                            drug_interaction_dataset,
                            cancer_types,
                            int(edge['from'][1:])
                        ))

                # get unique cancer genes out of all found
                ds_closest_baits = get_unique_bait_distances(ds_closest_baits)

                closest_baits = []
                minimum_distance = None

                # find closest bait out of baits
                if len(ds_closest_baits) > 0:
                    minimum_distance = min([cb['distance'] for cb in ds_closest_baits])
                    closest_baits = [cb for cb in ds_closest_baits if cb['distance'] <= minimum_distance]
                details_s['closest_cancer_genes'] = ','.join([cb['gene_name'] for cb in closest_baits])
                details_s['closest_distance'] = minimum_distance + 1 if minimum_distance else ''

            if node_graphId not in node_details:
                node_details[node_graphId] = details_s
            else:
                # TODO add node type
                node_details[node_graphId] = details_s

            if node_type == 'Node':
                genes.append(details_s)
            elif node_type == 'CancerNode':
                cancer_genes.append(details_s)
            elif node_type == 'Drug':
                drugs.append(details_s)

    # todo remove
    result['node_objects'] = genes
    result['cancerNode_objects'] = cancer_genes
    result['drug_objects'] = drugs

    # fetch information for used datasets
    result['cancer_dataset'] = serializers.CancerDatasetSerializer().to_representation(
        models.CancerDataset.objects.get(name=cancer_dataset))
    result['gene_interaction_dataset'] = serializers.InteractionDatasetSerializer().to_representation(
        models.InteractionGeneGeneDataset.objects.get(name=gene_interaction_dataset))
    result['drug_interaction_dataset'] = serializers.InteractionDatasetSerializer().to_representation(
        models.InteractionGeneDrugDataset.objects.get(name=drug_interaction_dataset))
    result['cancer_types'] = serializers.CancerTypeSerializer(many=True).to_representation(
        models.CancerType.objects.filter(id__in=cancer_types))

    if not view:
        return Response(result)
    else:
        if view == 'genes':
            if fmt == 'csv':
                items = []
                for i in genes:
                    new_i = {
                        'name': i['name'],
                        'entrez_id': i['entrez_id'],
                        'protein_name': i['protein_name'],
                        'uniprot_ac': i['uniprot_ac'],
                        'seed': is_seed[i['graphId']],
                        'interactome_degree': db_degrees[i['graphId']] if db_degrees else None
                        # 'closest_cancer_genes': i['closest_cancer_genes'],
                        # 'closest_distance': i['closest_distance']
                    }
                    if i.get('score'):
                        new_i['score'] = i['score']
                    items.append(new_i)
            else:
                items = genes
        elif view == 'cancer_driver_genes':
            if fmt == 'csv':
                items = []
                for i in cancer_genes:
                    new_i = {
                        'name': i['name'],
                        'entrez_id': i['entrez_id'],
                        'cancer_type': i['type'],
                        'protein_name': i['protein_name'],
                        'uniprot_ac': i['uniprot_ac'],
                        'seed': is_seed[i['graphId']],
                        'interactome_degree': db_degrees[i['graphId']] if db_degrees else None
                    }
                    if i.get('score'):
                        new_i['score'] = i['score']
                    items.append(new_i)
            else:
                items = cancer_genes
        elif view == 'drugs':
            if fmt == 'csv':
                items = []
                for i in drugs:
                    i['interactome_degree'] = db_degrees[i['graphId']] if db_degrees else None
                    if 'closest_distance' in i:
                        del i['closest_distance']
                    if 'closest_cancer_genes' in i:
                        del i['closest_cancer_genes']
                    items.append(i)
            else:
                items = drugs
        else:
            return Response({})
        if not fmt or fmt == 'json':
            return Response(items)
        elif fmt == 'csv':
            if len(items) != 0:
                keys = items[0].keys()
            else:
                keys = []
            response = HttpResponse(content_type='text/csv')
            response['Content-Disposition'] = 'attachment; filename='
            response['Content-Disposition'] += f'"{task.id}_{task.algorithm}_{gene_interaction_dataset}'
            response['Content-Disposition'] += f'_{drug_interaction_dataset}_{view}.csv"'
            dict_writer = csv.DictWriter(response, keys)
            dict_writer.writeheader()
            dict_writer.writerows(items)
            return response
        else:
            return Response({})


@api_view()
def graph_export(request):
    token_str = request.query_params['token']
    task = models.Task.objects.get(token=token_str)

    result = task_result(task)
    node_attributes = result.get('node_attributes')

    if not node_attributes:
        node_attributes = {}
        result['node_attributes'] = node_attributes

    network = result['network']
    cancer_types = result['cancer_types']
    cancer_types = list(map(int, cancer_types))

    node_types = node_attributes.get('node_types')
    if not node_types:
        node_types = {}
        node_attributes['node_types'] = node_types

    is_seed = node_attributes.get('is_seed')
    if not is_seed:
        is_seed = {}
        node_attributes['is_seed'] = is_seed

    scores = node_attributes.get('scores', {})

    node_details = {}
    node_attributes['details'] = node_details

    parameters = json.loads(task.parameters)
    seeds = parameters['seeds']
    nodes = network['nodes']
    edges = network['edges']

    G = nx.Graph()
    graph_id_to_name = {}
    for node in nodes:
        node_backend_id = node[1:]
        is_seed[node] = node in seeds
        node_type = node_types.get(node)
        details = None
        details_s = None
        if not node_type:
            node_type, details = infer_node_type_and_details(node)
        else:
            if node_type == 'Node':
                details = models.Gene.objects.get(id=node_backend_id)
            elif node_type == 'CancerNode':
                details = models.CancerGeneEntity.objects.filter(
                    cancer_gene_id=node_backend_id,
                    cancer_type_id__in=cancer_types
                ).first()
                # TODO: Add dataset to task and choose correct effect
            elif node_type == 'Drug':
                node_type = 'Drug'
                details = models.DrugEntity.objects.get(drug_id=node_backend_id)

        if not node_type or not details:
            continue

        if node_type == 'Node':
            details_s = serializers.GeneSerializer().to_representation(details)
        elif node_type == 'CancerNode':
            details_s = serializers.CancerGeneEntitySerializer().to_representation(details)
            # causes error in add_node function, rename key 'type'
            details_s['cancerType'] = details_s['type']
            del details_s['type']
        elif node_type == 'Drug':
            details_s = serializers.DrugEntitySerializer().to_representation(details)
            details_s.pop('trial_links', None)

        node_types[node] = node_type
        if scores.get(node) is not None:
            details_s['score'] = scores.get(node, None)
        node_details[node] = details_s
        if node_type == 'Node':
            G.add_node(details_s['name'], type=node_type, **details_s)
        elif node_type == 'CancerNode':
            G.add_node(details_s['name'], type=node_type, **details_s)
        elif node_type == 'Drug':
            G.add_node(details_s['name'], type=node_type, **details_s)
        graph_id_to_name[node] = details_s['name']

    for e in edges:
        G.add_edge(graph_id_to_name[e['from']], graph_id_to_name[e['to']])

    # replace missing attributes (e.g. some uniprot ACs that are missing) with an empty string, since
    # "write_graphml" does not support None
    for node in G.nodes:
        for attrib in G.nodes[node]:
            if G.nodes[node][attrib] is None:
                G.nodes[node][attrib] = ''

    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = f'attachment; filename="{task.id}_network.graphml"'
    nx.write_graphml(G, response)
    return response


class GeneDrugInteractionView(APIView):
    """
    Gene-Drug-Interaction Network
    """

    def get(self, request):
        if request.query_params.get('nodes'):
            backend_id_list = json.loads(request.query_params.get('nodes'))
            genes = list(models.Gene.objects.filter(id__in=backend_id_list).all())
        else:
            genes = []
            token = json.loads(request.query_params['token'])
            task = models.Task.objects.get(token=token)
            result = task_result(task)
            network = result['network']
            cancer_types = result['cancer_types']
            cancer_types = list(map(int, cancer_types))

            # only fetch drugs with selected status, -1 is all drugs
            drug_status = json.loads(request.query_params['drugStatus'])

            # cancer_dataset = result.get('cancer_dataset')
            # gene_interaction_dataset = result.get('gene_interaction_dataset')
            drug_interaction_dataset = result.get('drug_interaction_dataset')
            node_attributes = result.get('node_attributes')
            if not node_attributes:
                node_attributes = {}
            node_types = node_attributes.get('node_types')
            if not node_types:
                node_types = {}
            nodes = network['nodes']
            for node in nodes:
                node_type = node_types.get(node)
                details = None

                if not node_type:
                    node_type, details = infer_node_type_and_details(node, cancer_types)
                if node_type != 'Drug':  # should never occur since we searched for drug targets
                    if details is not None:
                        genes.append(details)
                    else:
                        try:
                            genes.append(models.Gene.objects.get(id=node[1:]))
                        except models.Gene.DoesNotExist:
                            pass

        pd_interactions = []
        drugs = set()
        genen_drug_interaction_object_list = models.GeneDrugInteraction.objects.filter(
            gene_id__in=genes,
            interaction_dataset_id__name=drug_interaction_dataset
        )
        for pdi_object in genen_drug_interaction_object_list:
            pd_interactions.append(pdi_object)
            drug = pdi_object.drug_id
            # filter for drug status
            if drug_status == -1 or drug.status.id == drug_status:
                drugs.add(drug)

        protein_drug_edges = {
            'genes': serializers.GeneSerializer(many=True).to_representation(genes),
            'drugs': serializers.DrugSerializer(many=True).to_representation(drugs),
            'edges': serializers.GeneDrugInteractionSerializer(many=True).to_representation(pd_interactions),
        }
        return Response(protein_drug_edges)


@api_view(['POST'])
def query_nodes(request):

    if request.data.get('nodes'):
        # nodes can be either gene entrezID, uniprot ID or drug db ID
        nodes = request.data.get('nodes')
    else:
        return HttpResponseBadRequest('nodes is missing')

    cancer_dataset_id = request.data.get('cancer_dataset', None)

    cancer_type_ids = request.data.get('cancer_types', '')
    # cancer_type_ids is a stringify-ed list of cancer type ids
    cancer_type_ids_list = cancer_type_ids.split(',')

    found_genes = []
    found_cancer_genes = []
    not_found = []
    nodes_genes = []
    not_known_cancer_genes = []  # alert user if genes are not known cancer genes
    for n in nodes:
        # if is drug db_id, find all related genes and append them for further search to input nodes
        try:
            drug_interactions = models.GeneDrugInteraction.objects.filter(drug__drug_id__db_id=n)
            for di in drug_interactions:
                nodes_genes.append(di.gene_id.entrez_id)
        except Exception:
            nodes_genes.append(n)

    for n in nodes_genes:
        if cancer_dataset_id is not None:
            # check if we have genes as cancer gene for given dataset

            # test if n is entrez ID
            # entrez ID will be integer (or float)
            cancer_gene = None
            if not str(n).isnumeric():
                n = str(n).upper()
                # test if node is uniprot ac
                result = models.CancerGeneEntity.objects.filter(
                    cancer_gene_id__uniprot_ac=n,
                    cancer_dataset_id=cancer_dataset_id,
                    cancer_type_id__in=cancer_type_ids_list
                )
                if result:
                    cancer_gene = result[0]
                else:
                    result = models.CancerGeneEntity.objects.filter(
                        cancer_gene_id__name=n,
                        cancer_dataset_id=cancer_dataset_id,
                        cancer_type_id__in=cancer_type_ids_list
                    )
                    if result:
                        cancer_gene = result[0]

            else:
                # node is not string, might be entrez ID
                result = models.CancerGeneEntity.objects.filter(
                    cancer_gene_id__entrez_id=int(n),
                    cancer_dataset_id=cancer_dataset_id,
                    cancer_type_id__in=cancer_type_ids_list
                )
                if result:
                    cancer_gene = result[0]

            if cancer_gene is not None:
                found_cancer_genes.append(serializers.CancerGeneEntitySerializer().to_representation(cancer_gene))
                continue

        # node has not been found in the cancer driver gene list, look in all genes
        # check if node is entrez id or uniprot_ac
        gene = None
        if not str(n).isnumeric():
            n = str(n).upper()
            result = models.Gene.objects.filter(
                uniprot_ac=n
            )
            if result:
                gene = result[0]
            else:
                result = models.Gene.objects.filter(
                    name=n
                )
                if result:
                    gene = result[0]
        else:
            # n is number, might be entrez
            result = models.Gene.objects.filter(
                entrez_id=int(n)
            )
            if result:
                gene = result[0]

        if gene is not None:
            gene_serialised = serializers.GeneSerializer().to_representation(gene)
            found_genes.append(gene_serialised)

            # check if genes are known cancer genes in other datasets
            cancer_gene_hits = models.CancerGeneEntity.objects.filter(cancer_gene_id=gene)
            if not len(cancer_gene_hits):
                not_known_cancer_genes.append(n)
            continue

        not_found.append(n)

    return Response({
        'genes': found_genes,
        'cancerGenes': found_cancer_genes,
        'notFound': not_found,
        'notKnownCancerGenes': not_known_cancer_genes
    })


class MutationScoreView(APIView):
    """
    mutation data of genes in mutation cancer types.
    """
    def post(self, request):
        if request.data.get('mutation_cancer_type', None) is None:
            return HttpResponseBadRequest('mutation_cancer_type is missing')

        mutation_cancer_type_object = models.MutationCancerType.objects.get(
            id=request.data.get('mutation_cancer_type')
        )

        if request.data.get('genes', ''):
            gene_ids = json.loads(request.data.get('genes', '')).split(',')
        else:
            return HttpResponseBadRequest('genes is missing')

        if request.data.get('cancer_genes', ''):
            cancer_gene_ids = json.loads(request.data.get('cancer_genes', '')).split(',')
        else:
            return HttpResponseBadRequest('cancer_genes is missing')

        cancer_driver_genes = []
        for gene_id in cancer_gene_ids:
            mutation_cancer_type_gene_object = models.MutationCounts.objects.filter(
                cancer_type=mutation_cancer_type_object,
                gene__id=gene_id
            ).first()

            if mutation_cancer_type_gene_object is None:
                # no information for this gene in this cancer type, set it to 0
                gene_object = models.Gene.objects.get(
                    id=gene_id
                )
                mutation_cancer_type_gene = serializers.GeneSerializer().to_representation(gene_object)
                mutation_cancer_type_gene['mutation_counts'] = None
                mutation_cancer_type_gene['mutation_score'] = None
            else:
                mutation_cancer_type_gene = serializers.MutationCountsGeneSerializer()\
                    .to_representation(mutation_cancer_type_gene_object)
            cancer_driver_genes.append(mutation_cancer_type_gene)
        genes = []
        for gene_id in gene_ids:
            mutation_cancer_type_gene_object = models.MutationCounts.objects.filter(
                cancer_type=mutation_cancer_type_object,
                gene__id=gene_id
            ).first()

            if mutation_cancer_type_gene_object is None:
                # no information for this gene in this cancer type, set it to 0
                gene_object = models.Gene.objects.get(
                    id=gene_id
                )
                mutation_cancer_type_gene = serializers.GeneSerializer().to_representation(gene_object)
                mutation_cancer_type_gene['mutation_counts'] = None
                mutation_cancer_type_gene['mutation_score'] = None
            else:
                mutation_cancer_type_gene = serializers.MutationCountsGeneSerializer() \
                    .to_representation(mutation_cancer_type_gene_object)

            genes.append(mutation_cancer_type_gene)
        return Response({'nodes': genes, 'cancerNodes': cancer_driver_genes})


class GeneExpressionView(APIView):
    """
    Expression of proteins in tissues.
    """

    def post(self, request):
        expression_cancer_type = models.ExpressionCancerType.objects.get(id=request.data.get('cancer_type'))

        if request.data.get('genes') and request.data.get('cancer_genes'):
            # filter genes and cancer genes separately in order to not have to sort them again in the frontend
            gene_backend_id_list = json.loads(request.data.get('genes', '')).split(',')
            if gene_backend_id_list[0] == '':
                genes = []
            else:
                genes = list(models.Gene.objects.filter(id__in=gene_backend_id_list).all())

            cancer_gene_backendId_list = json.loads(request.data.get('cancer_genes', '')).split(',')
            if cancer_gene_backendId_list[0] == '':
                cancer_genes = []
            else:
                cancer_genes = list(models.Gene.objects.filter(id__in=cancer_gene_backendId_list).all())

        elif request.data.get('data'):
            dataset = json.loads(request.data.get('data', '1'))
            cancer_types = json.loads(request.data.get('cancer_Types', '')).split(',')
            cancer_genes = []

            # cancer_dataset = models.CancerDataset.objects.get(id=dataset)
            cancer_genes.extend(
                models.CancerGeneEntity.objects.filter(
                    cancer_dataset_id=dataset, cancer_type_id__in=cancer_types
                     ).all()
            )

            cancer_gene_ids = [gene.id for gene in cancer_genes]

            genes = []
            edge_objects = models.GeneGeneInteraction.objects.filter(
                Q(gene_a__in=cancer_gene_ids) | Q(gene_b__in=cancer_gene_ids)
            )
            for edge_object in edge_objects:
                if edge_object.gene_a not in genes and edge_object.gene_a.id not in cancer_gene_ids:
                    genes.append(edge_object.gene_a)
                if edge_object.gene_b not in genes and edge_object.gene_b.id not in cancer_gene_ids:
                    genes.append(edge_object.gene_b)

        elif request.data.get('token'):
            genes = []
            cancer_genes = []
            task = models.Task.objects.get(token=request.data['token'])
            result = task_result(task)
            network = result['network']
            node_attributes = result.get('node_attributes')
            if not node_attributes:
                node_attributes = {}
            node_types = node_attributes.get('node_types')
            if not node_types:
                node_types = {}
            parameters = json.loads(task.parameters)
            seeds = parameters['seeds']
            nodes = network['nodes']
            for node in nodes + seeds:
                node_type = node_types.get(node)
                details = None
                if not node_type:
                    node_type, details = infer_node_type_and_details(node)
                if node_type == 'Node' or node_type == 'CancerNode':
                    if details:
                        genes.append(details)
                    else:
                        try:
                            gene = models.Gene.objects.get(uniprot_ac=node)
                            if gene not in genes:
                                genes.append(models.Gene.objects.get(uniprot_ac=node))
                        except models.Gene.DoesNotExist:
                            pass

        pt_expressions_genes = []
        for gene in genes:
            try:
                expression_level = models.GeneExpressionLevel.objects.get(
                    gene=gene,
                    expression_cancer_type=expression_cancer_type
                    )
                pt_expressions_genes.append({
                    'gene': serializers.GeneSerializer().to_representation(gene),
                    'level': expression_level.expression_level_tpm,
                    'score': expression_level.expression_level,
                })
            except models.GeneExpressionLevel.DoesNotExist:
                pt_expressions_genes.append({
                    'gene': serializers.GeneSerializer().to_representation(gene),
                    'level': None,
                    'score': None
                })

        pt_expressions_cancer_genes = []
        for gene in cancer_genes:
            try:
                expression_level = models.GeneExpressionLevel.objects.get(
                    gene=gene,
                    expression_cancer_type=expression_cancer_type
                    )
                pt_expressions_cancer_genes.append({
                    'gene': serializers.GeneSerializer().to_representation(gene),
                    'level': expression_level.expression_level_tpm,
                    'score': expression_level.expression_level,
                })
            except models.GeneExpressionLevel.DoesNotExist:
                pt_expressions_cancer_genes.append({
                    'gene': serializers.GeneSerializer().to_representation(gene),
                    'level': None,
                    'score': None
                })

        return Response({'genes': pt_expressions_genes, 'cancer_genes': pt_expressions_cancer_genes})


class TissueExpressionView(APIView):
    """
    Expression of proteins in tissues.
    """

    def post(self, request):
        tissue = models.Tissue.objects.get(id=request.data.get('tissue'))

        if request.data.get('genes') and request.data.get('cancer_genes'):
            # filter genes and cancer genes separately in order to not have to sort them again in the frontend
            gene_backend_id_list = json.loads(request.data.get('genes', '')).split(',')
            if gene_backend_id_list[0] == '':
                genes = []
            else:
                genes = list(models.Gene.objects.filter(id__in=gene_backend_id_list).all())

            cancer_gene_backendId_list = json.loads(request.data.get('cancer_genes', '')).split(',')
            if cancer_gene_backendId_list[0] == '':
                cancer_genes = []
            else:
                cancer_genes = list(models.Gene.objects.filter(id__in=cancer_gene_backendId_list).all())

        elif request.data.get('data'):
            dataset = json.loads(request.data.get('data', '1'))
            cancer_types = json.loads(request.data.get('cancer_types', '')).split(',')
            cancer_genes = []

            # cancer_dataset = models.CancerDataset.objects.get(id=dataset)
            cancer_genes.extend(
                models.CancerGeneEntity.objects.filter(
                    cancer_dataset_id=dataset, cancer_type_id__in=cancer_types
                     ).all()
            )

            cancer_gene_ids = [gene.id for gene in cancer_genes]

            genes = []
            edge_objects = models.GeneGeneInteraction.objects.filter(
                Q(gene_a__in=cancer_gene_ids) | Q(gene_b__in=cancer_gene_ids)
            )
            for edge_object in edge_objects:
                if edge_object.gene_a not in genes and edge_object.gene_a.id not in cancer_gene_ids:
                    genes.append(edge_object.gene_a)
                if edge_object.gene_b not in genes and edge_object.gene_b.id not in cancer_gene_ids:
                    genes.append(edge_object.gene_b)

        elif request.data.get('token'):
            genes = []
            cancer_genes = []
            task = models.Task.objects.get(token=request.data['token'])
            result = task_result(task)
            network = result['network']
            node_attributes = result.get('node_attributes')
            if not node_attributes:
                node_attributes = {}
            node_types = node_attributes.get('node_types')
            if not node_types:
                node_types = {}
            parameters = json.loads(task.parameters)
            seeds = parameters['seeds']
            nodes = network['nodes']
            for node in nodes + seeds:
                node_type = node_types.get(node)
                details = None
                if not node_type:
                    node_type, details = infer_node_type_and_details(node)
                if node_type == 'Node' or node_type == 'CancerNode':
                    if details:
                        genes.append(details)
                    else:
                        try:
                            gene = models.Gene.objects.get(uniprot_ac=node)
                            if gene not in genes:
                                genes.append(models.Gene.objects.get(uniprot_ac=node))
                        except models.Gene.DoesNotExist:
                            pass

        pt_expressions_genes = []
        for gene in genes:
            try:
                expression_level = models.ExpressionLevel.objects.get(gene=gene, tissue=tissue)
                pt_expressions_genes.append({
                    'gene': serializers.GeneSerializer().to_representation(gene),
                    'level': expression_level.expression_level,
                })
            except models.ExpressionLevel.DoesNotExist:
                pt_expressions_genes.append({
                    'gene': serializers.GeneSerializer().to_representation(gene),
                    'level': None,
                })

        pt_expressions_cancer_genes = []
        for gene in cancer_genes:
            try:
                expression_level = models.ExpressionLevel.objects.get(gene=gene, tissue=tissue)
                pt_expressions_cancer_genes.append({
                    'gene': serializers.GeneSerializer().to_representation(gene),
                    'level': expression_level.expression_level,
                })
            except models.ExpressionLevel.DoesNotExist:
                pt_expressions_cancer_genes.append({
                    'gene': serializers.GeneSerializer().to_representation(gene),
                    'level': None,
                })

        return Response({'genes': pt_expressions_genes, 'cancerGenes': pt_expressions_cancer_genes})


class TissueView(APIView):
    def get(self, request):
        tissues = models.Tissue.objects.all()
        return Response(serializers.TissueSerializer(many=True).to_representation(tissues))


class DrugTargetActionView(APIView):
    def get(self, request):
        # TODO use the model instead of strings
        # drug_target_actions = models.DrugTargetAction.objects.all()
        # return Response(
        #     serializers.DrugTargetActionSerializer(many=True).to_representation(drug_target_actions)
        #     )

        drug_target_actions = [
            {'name': 'inhibitor', 'backendId': 0},
            {'name': 'activator', 'backendId': 1},
            {'name': 'not inhibitor', 'backendId': 2},
            {'name': 'not activator', 'backendId': 3}
        ]
        return Response(drug_target_actions)


class ExpressionCancerTypeView(APIView):
    def get(self, request):
        expression_cancer_types = models.ExpressionCancerType.objects.all()
        return Response(
            serializers.ExpressionCancerTypeSerializer(many=True).to_representation(expression_cancer_types)
            )


class MutationCancerTypeView(APIView):
    def get(self, request):
        cancer_mutation_type_objects = models.MutationCancerType.objects.all()
        return Response(serializers.MutationCancerTypeSerializer(many=True)
                        .to_representation(cancer_mutation_type_objects))


class DrugStatusView(APIView):
    def get(self, request):
        drug_status = models.DrugStatus.objects.all()
        return Response(serializers.DrugStatusSerializer(many=True).to_representation(drug_status))


@api_view(['POST'])
def query_mutation_cancer_type_genes(request):
    threshold = request.data['threshold']
    mutation_cancer_type_id = request.data['mutation_cancer_type_id']
    cancer_types = request.data['cancer_types']
    mutation_cancer_type_object = models.MutationCancerType.objects.get(id=mutation_cancer_type_id)
    # get all genes that fit the threshold
    genes = []
    cts = models.MutationCounts.objects.filter(cancer_type=mutation_cancer_type_object, mutation_score__gte=threshold)
    for el in cts:
        genes.append(serializers.GeneSerializer().to_representation(el.gene))

    # find the cancer genes
    gene_ids = [gene['backendId'] for gene in genes]
    cancer_gene_objects = models.CancerGeneEntity.objects.filter(
        cancer_gene_id__in=gene_ids,
        cancer_type_id__in=cancer_types)
    cancer_genes = serializers.CancerGeneEntitySerializer(many=True).to_representation(cancer_gene_objects)

    # remove cancer genes from genes
    cancer_gene_ids = [gene['backendId'] for gene in cancer_genes]
    genes = [gene for gene in genes if gene['backendId'] not in cancer_gene_ids]

    return Response({
        'genes': genes,
        'cancerGenes': cancer_genes
    })


@api_view(['POST'])
def query_expression_cancer_type_genes(request):
    threshold = request.data['threshold']
    expression_cancer_type_id = request.data['expression_cancer_type_id']
    cancer_types = request.data['cancer_types']
    expression_cancer_type_object = models.ExpressionCancerType.objects.get(id=expression_cancer_type_id)
    # get all genes that fit the threshold
    genes = []
    expression_levels = models.GeneExpressionLevel.objects.filter(
        expression_cancer_type=expression_cancer_type_object,
        expression_level_tpm__gte=threshold)
    for el in expression_levels:
        genes.append(serializers.GeneSerializer().to_representation(el.gene))

    # find the cancer genes
    gene_ids = [gene['backendId'] for gene in genes]
    cancer_gene_objects = models.CancerGeneEntity.objects.filter(
        cancer_gene_id__in=gene_ids,
        cancer_type_id__in=cancer_types)
    cancer_genes = serializers.CancerGeneEntitySerializer(many=True).to_representation(cancer_gene_objects)

    # remove cancer genes from genes
    cancer_gene_ids = [gene['backendId'] for gene in cancer_genes]
    genes = [gene for gene in genes if gene['backendId'] not in cancer_gene_ids]

    return Response({
        'genes': genes,
        'cancerGenes': cancer_genes
    })


@api_view(['POST'])
def query_tissue_genes(request):
    threshold = request.data['threshold']
    tissue_id = request.data['tissue_id']
    cancer_types = request.data['cancer_types']
    tissue = models.Tissue.objects.get(id=tissue_id)
    # get all genes that fit the threshold
    genes = []
    for el in tissue.expressionlevel_set.filter(expression_level__gte=threshold):
        genes.append(serializers.GeneSerializer().to_representation(el.gene))

    # find the cancer genes
    gene_ids = [gene['backendId'] for gene in genes]
    cancer_gene_objects = models.CancerGeneEntity.objects.filter(
        cancer_gene_id__in=gene_ids,
        cancer_type_id__in=cancer_types)
    cancer_genes = serializers.CancerGeneEntitySerializer(many=True).to_representation(cancer_gene_objects)

    # remove cancer genes from genes
    cancer_gene_ids = [gene['backendId'] for gene in cancer_genes]
    genes = [gene for gene in genes if gene['backendId'] not in cancer_gene_ids]

    return Response({
        'genes': genes,
        'cancerGenes': cancer_genes
    })


class DiseaseView(APIView):
    def get(self, request):
        disease_objects = models.Disease.objects.all()
        diseases = serializers.DiseaseSerializer(many=True).to_representation(disease_objects)
        return Response({
            'diseases': diseases
        })


class GeneDrugInteractionLookup(APIView):
    def get(self, request):
        text = request.query_params.get('text')
        fmt = request.query_params.get('fmt') if 'fmt' in request.query_params else None
        dataset = request.query_params.get('dataset_id', None)

        if dataset is None:
            dataset_name = request.query_params.get('dataset')
            dataset = models.InteractionGeneDrugDataset.objects.get(name__iexact=dataset_name).id

        if text.isnumeric():
            gene_object = models.Gene.objects.filter(entrez_id=int(text)).first()
        else:
            gene_object = models.Gene.objects.filter(uniprot_ac=text).first()
            if gene_object is None:
                gene_object = models.Gene.objects.filter(name__iexact=text).first()
            if gene_object is None:
                gene_object = models.Gene.objects.filter(name__icontains=text).first()
            if gene_object is None:
                gene_object = models.Gene.objects.filter(protein_name__icontains=text).first()

        if gene_object is None:
            return Response({'found': False})

        # see if we have cancer gene, if so, fetch cancer types
        cancer_gene_objects = models.CancerGeneEntity.objects.filter(
            cancer_gene_id=gene_object
        )
        cancer_types = []
        cancer_datasets = {}
        for cancer_gene_object in cancer_gene_objects:
            cancer_type = serializers.CancerTypeSerializer().to_representation(cancer_gene_object.cancer_type_id)

            cancer_dataset = serializers.CancerDatasetSerializer().to_representation(
                cancer_gene_object.cancer_dataset_id)

            cancer_type['dataset'] = cancer_dataset['name']
            cancer_types.append(cancer_type)

            if cancer_gene_object.cancer_dataset_id not in cancer_datasets:
                cancer_datasets[cancer_gene_object.cancer_dataset_id] = cancer_dataset

        # get drug interactions of gene
        interactions = models.GeneDrugInteraction.objects.filter(
            gene_id=gene_object,
            interaction_dataset_id=dataset
        )
        drugs = [interaction.drug_id for interaction in interactions]
        drugs = serializers.DrugSerializer(many=True).to_representation(drugs)

        if fmt is None:
            gene = serializers.GeneSerializer().to_representation(gene_object)
            return Response({'found': True, 'gene': gene, 'drugs': drugs, 'cancer_types': cancer_types})
        elif fmt == 'csv':
            internal_keys = ['backendId', 'graphId']
            keys = list(drugs[0].keys()) if len(drugs) else []
            # remove internal data
            for k in internal_keys:
                keys.remove(k)

            items = []
            if len(keys):
                for item in drugs:
                    # remove internal data
                    for k in internal_keys:
                        del item[k]
                    items.append(item)

            response = HttpResponse(content_type='text/csv')
            response['Content-Disposition'] = f'attachment; filename="interaction_{gene_object.name}.csv"'
            dict_writer = csv.DictWriter(response, keys)
            dict_writer.writeheader()
            dict_writer.writerows(items)
            return response


class CancernetLookupView(APIView):
    def get(self, request):
        drug_id = request.query_params.get('drug_id')
        cancernet_objects = models.Cancernet.objects.filter(drug=int(drug_id))
        return Response({
            'cancernet': serializers.CancernetSerializer(many=True).to_representation(cancernet_objects)
        })


class DrugInteractionLookup(APIView):
    def get(self, request):
        text = request.query_params.get('text')
        fmt = request.query_params.get('fmt') if 'fmt' in request.query_params else None
        dataset = request.query_params.get('dataset_id', None)

        if dataset is None:
            dataset_name = request.query_params.get('dataset')
            dataset = models.InteractionGeneDrugDataset.objects.get(name__iexact=dataset_name).id

        if text[:2].upper() == 'DB':
            # is db id
            drug_object = models.Drug.objects.filter(db_id=text).first()
        else:
            # try to find name
            drug_object = models.Drug.objects.filter(name__iexact=text).first()

            if drug_object is None:
                drug_object = models.Drug.objects.filter(name__icontains=text).first()

        if drug_object is None:
            return Response({'found': False, 'genes': [], 'drug': []})

        interactions = models.GeneDrugInteraction.objects.filter(
            drug_id=drug_object,
            interaction_dataset_id=dataset
        )

        genes = []
        for interaction in interactions:
            gene = interaction.gene_id
            # see if we have cancer gene, if so, fetch cancer types
            cancer_gene_objects = models.CancerGeneEntity.objects.filter(
                cancer_gene_id=gene
            )

            cancer_types = []
            cancer_datasets = {}
            for cancer_gene_object in cancer_gene_objects:
                cancer_type = serializers.CancerTypeSerializer().to_representation(cancer_gene_object.cancer_type_id)

                cancer_dataset = serializers.CancerDatasetSerializer().to_representation(
                    cancer_gene_object.cancer_dataset_id)

                # add cancer gene dataset to cancer type name to make it more comprehensible
                cancer_type['name'] += f' ({cancer_dataset["name"]})'
                cancer_types.append(cancer_type)

                if cancer_gene_object.cancer_dataset_id not in cancer_datasets:
                    cancer_datasets[cancer_gene_object.cancer_dataset_id] = cancer_dataset

            genes.append({
                'gene': serializers.GeneSerializer().to_representation(gene),
                'cancerTypes': cancer_types,
                'cancerDatasets': cancer_datasets.values()
            })

        if fmt is None:
            drug = serializers.DrugSerializer().to_representation(drug_object)
            return Response({'found': True, 'genes': genes, 'drug': [drug]})
        elif fmt == 'csv':
            internal_keys = ['backendId', 'graphId', 'gene_alias']
            keys = list(genes[0]['gene'].keys()) if len(genes) else []
            # remove internal data
            for k in internal_keys:
                keys.remove(k)

            items = []
            if len(keys):
                keys.extend(['cancerTypes', 'cancerDatasets'])
                for gene in genes:
                    item = gene['gene']
                    item['cancerTypes'] = [x['name'] for x in gene['cancerTypes']]
                    item['cancerDatasets'] = [x['name'] for x in gene['cancerDatasets']]
                    # remove internal data
                    for k in internal_keys:
                        del item[k]

                    items.append(item)
            response = HttpResponse(content_type='text/csv')
            response['Content-Disposition'] = f'attachment; filename="interaction_{drug_object.name}.csv"'
            dict_writer = csv.DictWriter(response, keys)
            dict_writer.writeheader()
            dict_writer.writerows(items)

            return response


class DiseaseGeneInteractionView(APIView):

    def get(self, request):
        # get list of selected diseases
        diseases = json.loads(request.query_params.get('diseases'))

        # look up if genes are related to diseases
        if request.query_params.get('genes') and request.query_params.get('cancer_genes'):
            # nodes are given, look up which nodes occur in selected diseases
            gene_ids = json.loads(request.query_params.get('genes', []))
            cancer_gene_ids = json.loads(request.query_params.get('cancer_genes', []))

            # normal genes
            disease_gene_interactions = models.DiseaseGeneInteractions.objects.filter(
                gene__in=gene_ids,
                disease__in=diseases
            )
            in_diseases_genes = {gene_id: False for gene_id in gene_ids}
            for interaction in disease_gene_interactions:
                in_diseases_genes[f'g{interaction.gene.id}'] = True

            # cancer genes
            disease_cancer_gene_interactions = models.DiseaseGeneInteractions.objects.filter(
                gene__in=cancer_gene_ids,
                disease__in=diseases
            )
            in_diseases_cancer_genes = {gene_id: False for gene_id in cancer_gene_ids}
            for interaction in disease_cancer_gene_interactions:
                in_diseases_cancer_genes[f'g{interaction.gene.id}'] = True

            return Response({
                'in_diseases_genes': in_diseases_genes,
                'in_diseases_cancer_genes': in_diseases_cancer_genes
            })

        # fetch all genes in database for diseases
        elif request.query_params.get('data') and request.query_params.get('cancer_types'):
            cancer_dataset = json.loads(request.query_params.get('data'))
            cancer_types = json.loads(request.query_params.get('cancer_types'))
            # fetch all related genes for given diseases
            cancer_gene_objects = models.CancerGeneEntity.objects.filter(
                cancer_dataset_id=cancer_dataset,
                cancer_type_id__in=cancer_types
            )
            # lookup is id: cancer gene entity object
            cancer_genes_lookup = {g.cancer_gene_id.id: g for g in cancer_gene_objects}
            disease_gene_interactions = models.DiseaseGeneInteractions.objects.filter(
                disease__in=diseases
            )

            found_cancer_genes = []     # collect cancer gene entity objects
            found_genes = []    # collect gene ids, later we fetch gene objects
            for interaction in disease_gene_interactions:
                if interaction.gene.id in cancer_genes_lookup:
                    found_cancer_genes.append(cancer_genes_lookup[interaction.gene.id])
                else:
                    found_genes.append(interaction.gene.id)

            cancer_genes = serializers.CancerGeneEntitySerializer(many=True)\
                .to_representation(found_cancer_genes)

            gene_objects = models.Gene.objects.filter(id__in=found_genes)
            genes = serializers.GeneSerializer(many=True).to_representation(gene_objects)

            return Response({
                'genes': genes,
                'cancer_genes': cancer_genes
            })
