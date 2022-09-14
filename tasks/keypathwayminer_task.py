import base64
import json
import random
import string
import time
import datetime
from os.path import join

import requests

from caddie import models, serializers
from caddie.management.includes.DatabaseController import DatabaseController as DBC

from tasks.task_hook import TaskHook

# Base URL
url = 'http://172.25.0.1:9003/keypathwayminer/requests/'
# url = 'https://exbio.wzw.tum.de/keypathwayminer/requests/'
attached_to_id = ''.join(random.choices(string.ascii_uppercase + string.digits, k=32))


def send_request(sub_url, data):
    """
    Send a POST request with form-data to a given sub-URL and retrieve the JSON response.
    Throws a RuntimeError if there was an error while submitting

    :param sub_url: Sub-URL to send the POST request to
    :param data: Data dictionary that is sent via the POST request
    :return: JSON-object from the server request
    """
    request_url = join(url, sub_url)

    response = requests.post(url=request_url, data=data)


    # Check if submitting the job was successful
    if response.status_code != 200:
        raise RuntimeError(f'KPM server response code was "{response.status_code}", expected "200".')

    try:
        response_json = response.json()
    except json.decoder.JSONDecodeError:
        raise RuntimeError(f'The response could not be decoded as JSON, please check the URL:\n{request_url}')

    return response_json


def kpm_task(task_hook: TaskHook):
    """
    Run KeyPathwayMiner on given proteins and parameters remotely using the RESTful API of KPM-web
    Updates status of the TaskHook by polling the KPM-web server every second
    Writes results back to the TaskHook as 'networks'.

    :param task_hook: Needs to have 'k' set as a parameter (str or int) and a list of proteins set
    :return: None
    """
    cancer_dataset_id = task_hook.parameters.get("cancer_dataset", 1)
    cancer_types = cancer_types = task_hook.parameters.get("cancer_types", [])
    # --- Fetch and generate the datasets
    dataset_name = 'indicatorMatrix'
    indicator_matrix_string = ''

    # we need to convert the seeds from graphId into uniprot ac's
    seed_backend_ids = [seed_id[1:] for seed_id in task_hook.seeds]
    gene_objects = models.Gene.objects.filter(id__in=seed_backend_ids)
    gene_uniprot_acs = [gene.uniprot_ac for gene in gene_objects]

    for seed in gene_uniprot_acs:
        indicator_matrix_string += f'{seed}\t1\n'


    datasets = [
        {
            'name': dataset_name,
            'attachedToID': attached_to_id,
            'contentBase64': base64.b64encode(indicator_matrix_string.encode('UTF-8')).decode('ascii')
        }
    ]

    datasets_data = json.dumps(datasets)

    # --- Generate KPM settings
    k_val = str(task_hook.parameters['k'])
    kpm_settings = {
        'parameters': {
            'name': f'CADDIE run on {datetime.datetime.now()}',
            'algorithm': 'Greedy',
            'strategy': 'INES',
            'removeBENs': 'true',
            'unmapped_nodes': 'Add to negative list',
            'computed_pathways': 1,
            'graphID': 27,
            'l_samePercentage': 'false',
            'samePercentage_val': 0,
            'k_values': {
                'val': k_val,
                'val_step': '1',
                'val_max': k_val,
                'use_range': 'false',
                'isPercentage': 'false'
            },
            'l_values': {
                'val': '0',
                'val_step': '1',
                'val_max': '0',
                'use_range': 'false',
                'isPercentage': 'false',
                'datasetName': dataset_name
            }
        },
        'withPerturbation': 'false',
        'perturbation': [
            {
                'technique': 'Node-swap',
                'startPercent': '5',
                'stepPercent': '1',
                'maxPercent': '15',
                'graphsPerStep': '1'
            }
        ],
        'linkType': 'OR',
        'attachedToID': attached_to_id,
        'positiveNodes': '',
        'negativeNodes': ''
    }

    kpm_settings_data = json.dumps(kpm_settings)

    # --- Submit kpm job asynchronously
    kpm_job_data = {'kpmSettings': kpm_settings_data,
                    'datasets': datasets_data}

    submit_json = send_request('submitAsync', kpm_job_data)

    # Check if the submission was correct (add check whether parameters were correct)
    if not submit_json["success"]:
        raise RuntimeError(f'Job submission failed. Server response:\n{submit_json}')

    # Obtain questID for getting the result
    quest_id_data = {'questID': submit_json['questID']}
    # print(submit_json["resultUrl"])  # Remove in production

    # --- Retrieve status and update task_hook every 1s
    old_progress = -1
    while True:
        # Get status of job
        status_json = send_request('runStatus', quest_id_data)

        # Check if the questID exists (should)
        if not status_json['runExists']:
            raise RuntimeError(f'Job status retrieval failed. Run does not exist:\n{status_json}')

        # Set progress only when it changed
        progress = status_json['progress']
        if old_progress != progress:
            task_hook.set_progress(progress=progress, status='')
            old_progress = progress

        # Stop and go to results
        if status_json['completed'] or status_json['cancelled']:
            break

        time.sleep(1)

    # --- Retrieve results and write back
    results_json = send_request('results', quest_id_data)

    if not results_json['success']:
        raise RuntimeError(f'Job terminated but was unsuccessful:\n{results_json}')

    graphs_json = results_json.get('resultGraphs')

    # Build the networks
    network = None

    # Only build networks if the result is not empty
    if graphs_json is not None:

        # create 'is_seed' and 'node_types' for result loading
        is_seed = {}
        node_types = {}
        for graph in graphs_json:
            # Ignore the union set
            if graph['isUnionSet']:
                continue

            # Add nodes
            nodes = []
            for node in graph['nodes']:
                # we have the uniprotAc, fetch gene from database
                # see if we have cancer gene
                gene_object = models.CancerGeneEntity.objects.filter(
                    cancer_gene_id__uniprot_ac=node['name'],
                    cancer_dataset_id__name=cancer_dataset_id,
                    cancer_type_id__in=cancer_types
                ).first()
                gene_type = 'CancerNode'

                # if no cancer Gene is found, try to find normal gene
                if gene_object is None:
                    gene_object = DBC.get_gene(uniprot_ac=node['name'])
                    gene_type = 'Node'

                if gene_object is not None:
                    # node was found in db, add to result
                    if gene_type == 'CancerNode':
                        gene = serializers.CancerGeneEntitySerializer().to_representation(gene_object)
                    elif gene_type == 'Node':
                        gene = serializers.GeneSerializer().to_representation(gene_object)
                    nodes.append(gene['graphId'])

                    # add entry to 'node_types'
                    node_types[gene['graphId']] = gene_type
                    # add entry to 'is_seed'
                    if node['name'] in gene_uniprot_acs:
                        is_seed[gene['graphId']] = True
                    else:
                        is_seed[gene['graphId']] = False


            # Add edges
            edges = []
            for edge in graph['edges']:
                # get source and target genes
                source_gene_object = DBC.get_gene(uniprot_ac=edge['source'])
                target_gene_object = DBC.get_gene(uniprot_ac=edge['target'])

                if source_gene_object is not None and target_gene_object is not None:
                    # add edge to results if source and target gene are found in db
                    source_gene = serializers.GeneSerializer().to_representation(source_gene_object)
                    target_gene = serializers.GeneSerializer().to_representation(target_gene_object)

                    edges.append({'from': source_gene['graphId'], 'to': target_gene['graphId']})

            # Add nodes and edges to network
            network = {'nodes': nodes, 'edges': edges}

    result_dict = {
        'network': network,
        "node_attributes": {"node_types": node_types, "is_seed": is_seed, "db_degrees": {}, "is_result": {}},
        "cancer_types":  cancer_types,
        'cancer_dataset':  cancer_dataset_id,
        'gene_interaction_dataset': task_hook.parameters.get("gene_interaction_dataset", "BioGRID"),
        'drug_interaction_dataset': task_hook.parameters.get("drug_interaction_dataset", "BioGRID")
    }

    task_hook.set_results(results=result_dict)
