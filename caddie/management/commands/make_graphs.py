#!/usr/bin/env python
# coding: utf-8

import graph_tool.all as gt
from caddie import models
from caddie.management.includes.DatabaseController import DatabaseController
import multiprocessing
from django import db
from django.core.management import BaseCommand
import django

django.setup()

KERNEL = 1


def create_gt(params):
    """
    create a graph tool file with all genes and cancer genes for each dataset
    run this function only in docker since we have graph_tools installed in docker
    """
    dataset = params
    # get data from api
    print(f'genes/{dataset}')
    data = DatabaseController.internal_genes(dataset)

    print(f'edges')
    data['edges'] = DatabaseController.internal_gene_interactions_all()

    print(f'drugs')
    drug_data = DatabaseController.internal_drugs_all()
    data['drugs'] = drug_data['drugs']
    data['drugEdges'] = drug_data['drugEdges']

    # preprocessing: collect all cancer types for each cancer node
    # entrezId are key
    cancerNodes = {}
    for node in data['cancerNodes']:
        if node['graphId'] not in cancerNodes:
            t = node['type_id']
            node['type'] = [str(t)]
            cancerNodes[node['graphId']] = node
        else:
            cancerNodes[node['graphId']]['type'].append(str(node['type_id']))

    for node in cancerNodes:
        cancerNodes[node]['type'] = ','.join(cancerNodes[node]['type'])

    g = gt.Graph(directed=False)
    e_type = g.new_edge_property("string")
    e_backendId = g.new_edge_property("string")
    e_graphId = g.new_edge_property("string")
    e_cancer = g.new_edge_property("string")
    e_action = g.new_edge_property("string")
    e_dataset_name = g.new_edge_property("string")
    e_dataset_name_internal = g.new_edge_property("string")

    v_type = g.new_vertex_property("string")
    v_name = g.new_vertex_property("string")
    v_backendId = g.new_vertex_property("string")
    v_graphId = g.new_vertex_property("string")
    v_entrezId = g.new_vertex_property("string")
    v_cancer = g.new_vertex_property("string")
    # for drugs
    v_status = g.new_vertex_property("string")
    v_in_trial = g.new_vertex_property("string")
    v_in_literature = g.new_vertex_property("string")
    v_ctrpv2_id = g.new_vertex_property("int")
    v_db_id = g.new_vertex_property("string")
    v_mutation_scores = g.new_vertex_property("object")
    v_expression_scores = g.new_vertex_property("object")
    v_cancer_expression_scores = g.new_vertex_property("object")
    v_is_nutraceutical = g.new_vertex_property("bool")
    v_is_antineoplastic_and_immunomodulating_agent = g.new_vertex_property("bool")

    g.edge_properties["action"] = e_action
    g.edge_properties["type"] = e_type
    g.edge_properties["cancer"] = e_cancer
    g.edge_properties["dataset_name"] = e_dataset_name
    g.edge_properties["dataset_name_internal"] = e_dataset_name_internal
    g.edge_properties["backendId"] = e_backendId
    g.edge_properties["graphId"] = e_graphId

    g.vertex_properties["type"] = v_type
    g.vertex_properties["name"] = v_name
    g.vertex_properties["backendId"] = v_backendId
    g.vertex_properties["graphId"] = v_graphId
    g.vertex_properties["entrezId"] = v_entrezId
    g.vertex_properties["cancer"] = v_cancer
    g.vertex_properties["status"] = v_status
    g.vertex_properties["in_trial"] = v_in_trial
    g.vertex_properties["in_literature"] = v_in_literature
    g.vertex_properties["ctrpv2_id"] = v_ctrpv2_id
    g.vertex_properties["is_nutraceutical"] = v_is_nutraceutical
    g.vertex_properties["is_antineoplastic_and_immunomodulating_agent"] = v_is_antineoplastic_and_immunomodulating_agent
    g.vertex_properties["db_id"] = v_db_id
    g.vertex_properties["mutation_scores"] = v_mutation_scores
    g.vertex_properties["expression_scores"] = v_expression_scores
    g.vertex_properties["cancer_expression_scores"] = v_cancer_expression_scores

    # store nodes to connect them when creating edges
    vertices = {}
    drug_vertices = {}
    # add vertices
    print("adding nodes")
    # extend node data by cancer nodes, we create a normal node for each cancer node.
    # on reading the data, we decide which one to keep based on the user selected cancer types
    data['nodes'].extend(cancerNodes.values())
    for node in data['nodes']:
        v = g.add_vertex()
        v_type[v] = 'Node'
        v_name[v] = node['name']
        v_backendId[v] = node['backendId']
        v_graphId[v] = node['graphId']
        v_entrezId[v] = node['entrez_id']

        if node['graphId'] not in vertices:
            vertices[node['graphId']] = [int(v)]
        else:
            vertices[node['graphId']].append(int(v))

        # get scores for node
        scores = DatabaseController.internal_gene_scores(node['backendId'])
        v_mutation_scores[v] = scores['mutation_scores']
        v_expression_scores[v] = scores['expression_scores']
        v_cancer_expression_scores[v] = scores['cancer_expression_scores']

    print("done with nodes")

    print("adding cancer nodes")
    for backenId, node in cancerNodes.items():
        v = g.add_vertex()
        v_type[v] = 'CancerNode'
        v_name[v] = node['name']
        v_backendId[v] = node['backendId']
        v_graphId[v] = node['graphId']
        v_entrezId[v] = node['entrez_id']
        v_cancer[v] = node['type']

        if node['graphId'] not in vertices:
            vertices[node['graphId']] = [int(v)]
        else:
            vertices[node['graphId']].append(int(v))

        # get scores for node
        scores = DatabaseController.internal_gene_scores(node['backendId'])
        v_mutation_scores[v] = scores['mutation_scores']
        v_expression_scores[v] = scores['expression_scores']
        v_cancer_expression_scores[v] = scores['cancer_expression_scores']
    print("done with cancer nodes")

    print("adding drugs")
    for node in data['drugs']:
        v = g.add_vertex()
        v_type[v] = 'Drug'
        v_name[v] = node['name']
        v_status[v] = node['status']
        v_in_trial[v] = node['in_trial']
        v_in_literature[v] = node['in_literature']
        v_ctrpv2_id[v] = node['ctrpv2_id']
        v_backendId[v] = node['backendId']
        v_graphId[v] = node['graphId']
        v_db_id[v] = node['db_id']
        v_is_nutraceutical[v] = node['is_nutraceutical']
        v_is_antineoplastic_and_immunomodulating_agent[v] = node['is_atc_antineoplastic_and_immunomodulating_agent']
        drug_vertices[node['graphId']] = int(v)
    print("done with drugs")

    # add edges
    print("adding edges")
    for edge in data['edges']:
        a_indices = vertices[edge['interactor_a_graphId']]
        b_indices = vertices[edge['interactor_b_graphId']]

        done = set()
        for a in a_indices:
            for b in b_indices:
                if (a, b) and (b, a) not in done:
                    e = g.add_edge(a, b)
                    e_action[e] = ''
                    e_type[e] = 'gene-gene'
                    e_backendId[e] = edge['backendId']
                    e_graphId[e] = edge['graphId']
                    e_dataset_name[e] = edge['dataset_name']
                    e_dataset_name_internal[e] = edge['dataset_name_internal']

                    # get related cancer types
                    cancer_types = []
                    if edge['interactor_a_graphId'] in cancerNodes:
                        cancer_types.append(cancerNodes[edge['interactor_a_graphId']]['type'])
                    if edge['interactor_b_graphId'] in cancerNodes:
                        cancer_types.append(cancerNodes[edge['interactor_b_graphId']]['type'])
                    # sort out duplicates
                    # ['1,2,1', '1', '2,3']
                    s = ','.join(cancer_types)
                    splitted_string = s.split(',')
                    splitted_string = set(splitted_string)
                    e_cancer[e] = ','.join(splitted_string)

                    done.add((a, b))
                    done.add((b, a))
    print("done with edges")

    print("adding drug edges")
    for edge in data['drugEdges']:
        genes, drug = vertices[edge['gene_graphId']], drug_vertices[edge['drug_graphId']]
        for gene in genes:
            e = g.add_edge(drug, gene)
            e_action[e] = edge['action'] if edge['action'] else ''
            e_type[e] = 'drug-gene'
            e_backendId[e] = edge['backendId']
            e_graphId[e] = edge['graphId']
            e_dataset_name[e] = edge['dataset_name']
            e_dataset_name_internal[e] = edge['dataset_name_internal']

            # get cancer types edge is related to
            cancer_types = []
            if edge['gene_graphId'] in cancerNodes:
                cancer_types.append(cancerNodes[edge['gene_graphId']]['type'])
            e_cancer[e] = ','.join(cancer_types)

    print("done with drug edges")

    # save graph
    g.save(f"./data/networks/internal_{dataset}.gt")
    print(f"saved {dataset}")
    return g, cancerNodes, vertices


class Command(BaseCommand):
    def add_arguments(self, parser):
        pass

    def handle(self, *args, **kwargs):
        cancer_datasets = models.CancerDataset.objects.all()
        cancer_dataset_names = [e.name for e in cancer_datasets]

        parameter_combinations = []
        for cancer_dataset_name in cancer_dataset_names:
            parameter_combinations.append(cancer_dataset_name)

        # close all database connections so subprocesses will create their own connections
        # this prevents the processes from running into problems because of using the same connection
        db.connections.close_all()
        pool = multiprocessing.Pool(KERNEL)
        pool.map(create_gt, parameter_combinations)
        print('Made all graphs!')
        return
