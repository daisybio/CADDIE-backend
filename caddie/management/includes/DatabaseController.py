import logging
from django.db.models import Q
from django.db import OperationalError, IntegrityError
import caddie.models as models
import caddie.serializers as serializers
from django.db import connection
import pandas as pd
import numpy as np

logging.basicConfig(filename='logging.txt',
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)


class DatabaseController:
    """
    Interface to controll the postgresql database connected to django
    Provides a variety of functions to add and get data

    Is used to add data while populating in the beginning and also on runtime to retrieve data and #
    modify the database.
    Functions are called in the beginning by the DatabasePopulator.
    """

    def __init__(self):
        pass

    def delete_model(self, model):
        cursor = connection.cursor()
        try:
            cursor.execute(f'TRUNCATE TABLE "{model._meta.db_table}" RESTART IDENTITY CASCADE')
        except OperationalError:
            cursor.execute(f'DELETE FROM "{model._meta.db_table}"')
        # reset id autoincrementation
        # this is dangerous but needed due to the current handling of the graph_tool files
        cursor.execute(f'ALTER SEQUENCE "{model._meta.db_table}_id_seq" RESTART WITH 1')

    def delete_models(self, model_list):
        """
        Clears all models in the database listed in model_list.

        :param model_list: list(str) List of model names to delete
        :return:
        """
        for model_name in model_list:
            logging.info(f'Deleting {model_name} model ...')
            if model_name == 'DrugDataset':
                self.delete_model(models.DrugDataset)
            elif model_name == 'CancerDataset':
                self.delete_model(models.CancerDataset)
            elif model_name == 'CancerType':
                self.delete_model(models.CancerType)
            elif model_name == 'Drug':
                self.delete_model(models.Drug)
            elif model_name == 'DrugEntity':
                self.delete_model(models.DrugEntity)
            elif model_name == 'Gene':
                self.delete_model(models.Gene)
            elif model_name == 'RelationDatasetType':
                self.delete_model(models.RelationDatasetType)
            elif model_name == 'CancerGeneEntity':
                self.delete_model(models.CancerGeneEntity)
            elif model_name == 'GeneGeneInteraction':
                self.delete_model(models.GeneGeneInteraction)
            elif model_name == 'RelationGeneDrug':
                self.delete_model(models.GeneDrugInteraction)
            elif model_name == 'InteractionGeneGeneDataset':
                self.delete_model(models.InteractionGeneGeneDataset)
            elif model_name == 'InteractionGeneDrugDataset':
                self.delete_model(models.InteractionGeneDrugDataset)
            elif model_name == 'GeneDrugInteraction':
                self.delete_model(models.GeneDrugInteraction)
            elif model_name == 'Disease':
                self.delete_model(models.Disease)
            elif model_name == 'InteractionGeneDrugDataset':
                self.delete_model(models.InteractionGeneDrugDataset)
            elif model_name == 'DiseaseInteractions':
                self.delete_model(models.DiseaseInteractions)
            elif model_name == 'DiseaseGeneInteractions':
                self.delete_model(models.DiseaseGeneInteractions)
            elif model_name == 'ShortestDistanceGeneToCancerGene':
                self.delete_model(models.ShortestDistanceGeneToCancerGene)
            elif model_name == 'ShortestDistanceDrugToCancerGene':
                self.delete_model(models.ShortestDistanceDrugToCancerGene)
            elif model_name == 'MinSpanningTree':
                self.delete_model(models.MinSpanningTree)
            elif model_name == 'MutationCounts':
                self.delete_model(models.MutationCounts)
            elif model_name == 'MutationCancerType':
                self.delete_model(models.MutationCancerType)
            elif model_name == 'ExpressionCancerType':
                self.delete_model(models.ExpressionCancerType)
            elif model_name == 'GeneExpressionLevel':
                self.delete_model(models.GeneExpressionLevel)
            elif model_name == 'Tissue':
                self.delete_model(models.Tissue)
            elif model_name == 'ExpressionLevel':
                self.delete_model(models.ExpressionLevel)
            elif model_name == 'Cancernet':
                self.delete_model(models.Cancernet)

    def __count_gene_relations(self, gene_object):
        """
        Counts all gene-gene interactions for given gene in the database
        :param gene_object:
        :return:
        """
        # interactions with other genes
        gene_interactions = models.GeneGeneInteraction.objects.filter(Q(gene_a=gene_object) | Q(gene_b=gene_object))
        return len(gene_interactions)

    def add_cancer_type_occurrences(self):
        """
        Adds cancer_occurrences to the gene objects, information about many times gene occurs in cancer
        :return:
        """

        # iterate over all cancer datasets, information is just for each dataset and not globally
        cancer_datasets = models.CancerDataset.objects.all()
        for cancer_dataset in cancer_datasets:
            dataset_backendId = cancer_dataset.id

            # take all genes for dataset into account
            cancer_gene_objects = models.CancerGeneEntity.objects.filter(cancer_dataset_id=dataset_backendId)

            for cancer_gene_object in cancer_gene_objects:

                # 'cancer_dataset_id', 'cancer_type_id', 'cancer_gene_id' are unique together
                occurrences = models.CancerGeneEntity.objects\
                    .filter(cancer_dataset_id=dataset_backendId,
                            cancer_gene_id__id=cancer_gene_object.cancer_gene_id.id)\
                    .count()

                # set occurrences to cancer gene object
                cancer_gene_object.cancer_occurrences = occurrences
                cancer_gene_object.save()

    def __get_entrez_ids(self):
        """
        Fetches all entrez ids stored in Database
        :return: Gene[]
        """
        cancer_gene_objects = models.Gene.objects.all()

        cancer_genes = serializers.GeneSerializer(many=True).to_representation(cancer_gene_objects)

        cancer_genes_entrez_ids = [gene['entrez_id'] for gene in cancer_genes]

        return cancer_genes_entrez_ids

    def __get_cancer_type(self, cancer_type: str):
        """
        Returns the cancer_type id given a cancer_type string
        :param cancer_type:
        :return: cancer_type_object if in db else None
        """

        cancer_type_object_list = models.CancerType.objects.filter(name__iexact=cancer_type)

        if not len(cancer_type_object_list):
            # no cancer type was found for given name
            logging.debug('No cancer-type could be found for {}\n'.format(cancer_type))
            return None

        # list contains one item since search is unambiguously
        cancer_type_object = cancer_type_object_list[0]

        return cancer_type_object

    def __get_drug_dataset(self, drug_dataset: str):
        """
        Returns the drug_dataset id given a cancer_dataset string
        :param drug_dataset:
        :return: drug dataset object if in db else None
        """

        # filter case insensitive (name__iexact)
        drug_dataset_object_list = models.DrugDataset.objects.filter(name__iexact=drug_dataset)

        if not len(drug_dataset_object_list):
            # no drug dataset was found for given name
            logging.debug('No drug-dataset could be found for {}\n'.format(drug_dataset))
            return None

        # list contains one item since search is unambiguously
        drug_dataset_object = drug_dataset_object_list[0]

        return drug_dataset_object

    def __get_cancer_dataset(self, cancer_dataset: str):
        """
        Returns the cancer_dataset id given a cancer_dataset string
        :param cancer_dataset:
        :return: cancer dataset object if in db else None
        """

        cancer_dataset_object_list = models.CancerDataset.objects.filter(name__iexact=cancer_dataset)

        if not len(cancer_dataset_object_list):
            # no cancer dataset was found for given name
            logging.debug('No cancer-dataset could be found for {}\n'.format(cancer_dataset))
            return None

        # list contains one item since search is unambiguously
        cancer_dataset_object = cancer_dataset_object_list[0]

        return cancer_dataset_object

    def __get_disease(self, mondo_id: int = None, icd_10: str = None):
        """
        Returns the disease object for given mondo or icd_10 id
        :param mondo_id: MONDO-id like MONDO:1234567 (len = 13)
        :param icd10: ICD10 ID like A03
        :return: models.Disease
        """

        if mondo_id:
            try:
                return models.Disease.objects.filter(mondo_id=mondo_id)[0]
            # TODO find proper exception
            except Exception as e:
                logging.debug(f'Mondo ID {mondo_id} not found. {e}')
                return None

        if icd_10:
            try:
                return models.Disease.objects.filter(icd_10=icd_10)[0]
            # TODO find proper exception
            except Exception as e:
                logging.debug(f'ICD10 {icd_10} not found. {e}')
                return None

    @staticmethod
    def get_gene(entrez_id: int = None, uniprot_ac: str = None):
        """
        Returns the protein id given a cancer_dataset string
        :param entrez_id: int optional, either entrez_id or protein has to be passed
        :param gene: STRING optional, either entrez_id or protein has to be passed
        :return: gene_object if in db else None
        """
        gene_object = None

        if not entrez_id and not uniprot_ac:
            raise ValueError('ERROR: Either entrez_id or uniprot_ac has to be passed')

        if entrez_id is not None:

            gene_object = models.Gene.objects.filter(entrez_id=entrez_id).first()

            if gene_object is None:
                # no cancer protein was found
                logging.debug('No cancer-gene could be found for entrez_id: {}\n'.format(entrez_id))

                # TODO create own function for alias

        # uniprot has to be set here, but we still test it since runtime does not matter and code is more readable
        elif uniprot_ac is not None:
            gene_object = models.Gene.objects.filter(uniprot_ac=uniprot_ac).first()

            if gene_object is None:
                # no cancer gene was found
                logging.debug('No cancer-gene could be found for uniprot_ac: {}\n'.format(uniprot_ac))

        if gene_object is None:
            return None

        return gene_object

    def __get_cancer_gene_entity(self, cancer_dataset_id: int, cancer_type_id: int, cancer_gene_id: int):
        """
        Returns the entity_gene_id for given parameters
        :param cancer_dataset_id:
        :param cancer_type_id:
        :param cancer_gene_id:
        :return: id or None if non-existent
        """

        # attributes are unique, use get to retrieve exactly one row
        cancer_gene_entity_object_list = models.CancerGeneEntity.objects.filter(
            cancer_dataset_id=cancer_dataset_id,
            cancer_type_id=cancer_type_id,
            cancer_gene_id=cancer_gene_id
        )

        if not len(cancer_gene_entity_object_list):
            # no entity protein was found
            logging.debug('No cancer-protein could be found for cancer_dataset_id: {}, cancer_type_id: {}, '
                          'cancer_gene_id: {}\n'.format(cancer_dataset_id, cancer_type_id, cancer_gene_id))
            return None

        # list contains one item since search is unambiguously
        cancer_gene_entity_object = cancer_gene_entity_object_list[0]

        return cancer_gene_entity_object

    def __get_gene_gene_relation(self, gene_object_a: models.Gene, gene_object_b: models.Gene,
                                 interaction_dataset_object: models.GeneGeneInteraction):
        """
        returns a protein protein relation object if it exists, else None
        :param entrez_a:
        :param entrez_b:
        :return: gene_gene_relation_object or None
        """
        gene_gene_relation_object_list = models.GeneGeneInteraction.objects.filter(
            gene_a=gene_object_a,
            gene_b=gene_object_b,
            interaction_dataset_id=interaction_dataset_object
        )

        if not len(gene_gene_relation_object_list):
            # no interaction was found
            logging.debug('No interaction could be found for genes {} and {}\n'.format(
                gene_object_a.id, gene_object_b.id))
            return None

        # list contains one item since search is unambiguously
        gene_gene_relation_object = gene_gene_relation_object_list[0]

        return gene_gene_relation_object

    def __get_drug(self, name: str = None, db_id: str = None):
        """
        Returns the drug object given a drug name, if it does not exist in db returns None
        :param db_id:
        :param name:
        :return:
        """
        if db_id:
            drug_object_list = models.Drug.objects.filter(db_id=db_id)
        elif name:
            # utf8 encoding
            name = name.encode()
            drug_object_list = models.Drug.objects.filter(name=name)

        if not len(drug_object_list):
            # no drug was found
            logging.debug('No drug could be found for drug-name {} or db_id {}\n'.format(name, db_id))
            return None

        # list contains one item since search is unambiguously
        drug_object = drug_object_list[0]

        return drug_object

    def __get_drug_gene_interaction(self, drug, gene, dataset):
        """
        Returns the drug protein interaction object given a drug name, if it does not exist in db returns None
        :param drug: Drug
        :param gene: Gene
        :return:
        """

        drug_gene_interaction_object_list = models.GeneDrugInteraction.objects.filter(
            gene_id=gene.id,
            drug_id=drug.id,
            interaction_dataset_id=dataset
        )

        if not len(drug_gene_interaction_object_list):
            # no drug was found
            logging.debug('No drug protein interaction could be found for protein {} and drug {}\n'.format(
                gene.id, drug.id))
            return None

        # list contains one item since search is unambiguously
        drug_gene_interaction_object = drug_gene_interaction_object_list[0]

        return drug_gene_interaction_object

    def __add_cancer_type_dataset_relation(self, cancer_type: models.CancerType, cancer_dataset: models.CancerDataset):
        """
        Add a cancer type dataset relation
        :param cancer_type: string
        :return: cancer_type_id. id of cancer_type that got created or None if cancer-type already exists
        """
        try:
            cancer_type_dataset_relation_object = models.RelationDatasetType.objects.create(
                cancer_dataset_id=cancer_dataset,
                cancer_type_id=cancer_type
            )

            logging.debug('Created entry for cancer type dataset relation {}\n'.format(cancer_type))
            return cancer_type_dataset_relation_object
        except IntegrityError as e:
            # cancer type already exists
            logging.debug(e)
            logging.debug('Cancer cancer type dataset relatio {} might already exist (see message above)\n'
                          .format(cancer_type))
            return None

    def __add_cancer_type(self, cancer_type: str):
        """
        Add a cancer type to the database
        :param cancer_type: string
        :return: cancer_type_id. id of cancer_type that got created or None if cancer-type already exists
        """
        if not isinstance(cancer_type, str):
            raise TypeError('cancer_type must be a string')

        try:
            cancer_type_object = models.CancerType.objects.create(name=cancer_type)

            logging.debug('Created entry for cancer type {}\n'.format(cancer_type))
            return cancer_type_object
        except IntegrityError as e:
            # cancer type already exists
            logging.debug(e)
            logging.debug('Cancer type {} might already exist (see message above)\n'.format(cancer_type))
            return None

    def __add_cancer_dataset(self, dataset_name, link, version):
        """
        Adds a cancer_dataset entry if name is not in table
        :param dataset_name:
        :return: cancer_dataset or None if entry already exists
        """
        if not isinstance(dataset_name, str):
            raise TypeError('cancer_dataset must be a string')

        try:
            cancer_dataset_object = models.CancerDataset.objects.create(name=dataset_name, link=link, version=version)

            logging.debug('Created entry for cancer dataset {}\n'.format(dataset_name))
            return cancer_dataset_object

        except IntegrityError as e:
            # cancer type already exists
            logging.debug(e)
            logging.debug('Cancer dataset {} might already exist (see message above)\n'.format(dataset_name))
            return None

    def __add_gene(self, gene, entrez_id=None, alias='', protein_name=None, uniprot_ac=None):
        """
        Adding a protein and a cancer_type to the database.
        Creates entry in "genes" table
        :param entrez_id: INT entrez_id of protein
        :param gene: STRING name of cancer protein MANDATORY
        :param alias: gene alias, comma separated string
        :param protein_name: string for protein name
        :param uniprot_ac: string for uniprot ac
        :return: cancer_gene_object (Gene) or NONE if gene already in table
        """

        try:
            if pd.isnull(entrez_id):
                entrez_id = None

            cancer_gene_object = models.Gene.objects.create(
                name=gene,
                entrez_id=entrez_id,
                alias=alias,
                protein_name=protein_name,
                uniprot_ac=uniprot_ac
            )

            logging.debug('Created entry for gene {}\n'.format(gene))
            return cancer_gene_object

        except IntegrityError as e:
            # cancer type already exists
            logging.debug(e)
            logging.debug('Cancer gene {} might already exist (see message above)\n'.format(gene))
            return None

    def __add_cancer_gene_entity_instance(self, cancer_dataset_object, cancer_type_object, cancer_gene_object,
                                          pubmed_id=None):
        """
        Adds an entry to the table "cancer_gene_entities"
        :param pubmed_id: INT or NONE
        :param cancer_gene_object: INT
        :param cancer_type_object: INT
        :param cancer_dataset_object: INT
        :return: entity_id if entity-entry added to table or False if something goes wrong
        (entry probably already exists)
        """

        assert isinstance(cancer_gene_object, models.Gene)
        assert isinstance(cancer_type_object, models.CancerType)
        assert isinstance(cancer_dataset_object, models.CancerDataset)
        assert isinstance(pubmed_id, int) or pubmed_id is None

        try:
            cancer_gene_entity_object = models.CancerGeneEntity.objects.create(
                cancer_dataset_id=cancer_dataset_object,
                cancer_type_id=cancer_type_object,
                cancer_gene_id=cancer_gene_object,
                pubmed_id=pubmed_id
            )

            logging.debug('Created entry for gene entity id: {}, cancer_dataset_id: {}, '
                          'cancer_type_id: {}, cancer_gene_id: {}\n'
                          .format(cancer_gene_entity_object.id,
                                  cancer_dataset_object.id,
                                  cancer_type_object.id,
                                  cancer_gene_object.id))
            return cancer_gene_entity_object

        except IntegrityError as e:
            # cancer_gene_entity_object already exists
            logging.debug(e)
            logging.debug('ERROR: Entry for gene entity already exists, cancer_dataset_id: {}, '
                          'cancer_type_id: {}, cancer_gene_id: {}\n'
                          .format(cancer_dataset_object.id, cancer_type_object.id, cancer_gene_object.id))
            return None

    def __add_gene_gene_relation(self, gene_a, gene_b, interaction_dataset_object, data=None):
        """
        :param data: dict-like object with data
        :return:
        """

        if data is None:
            data = {
                'pubmed_id': '',
                'confidence_value': '',
                'detection_method_psi_mi': '',
                'detection_method_name': '',
                'type_psi_mi': '',
                'type_name': ''
            }

        try:

            if 'pubmed_id' in data:
                gene_gene_interaction_object = models.GeneGeneInteraction.objects.create(
                    gene_a=gene_a,
                    gene_b=gene_b,
                    interaction_dataset_id=interaction_dataset_object,
                    pubmed_id=data['pubmed_id'],
                    confidence_value=data['confidence_value'],
                    detection_method_psi_mi=data['detection_method_psi_mi'],
                    detection_method_name=data['detection_method_name'],
                    type_psi_mi=data['type_psi_mi'],
                    type_name=data['type_name']
                )
            else:
                gene_gene_interaction_object = models.GeneGeneInteraction.objects.create(
                    gene_a=gene_a,
                    gene_b=gene_b,
                    interaction_dataset_id=interaction_dataset_object
                )

            logging.debug('Created gene gene relation for genes {} and {} with id {}\n'.format(
                gene_a.id,
                gene_b.id,
                gene_gene_interaction_object.id
            ))
            return gene_gene_interaction_object

        except IntegrityError as e:
            # cancer type probably exists
            logging.debug(e)
            logging.debug(
                'Could not create gene gene relation for genes {} and {}\n'.format(gene_a.id, gene_b.id))
            return

    def __add_drug_gene_relation(self, drug, gene, interaction_dataset_object, pubmed_id=None):
        """
        :param gene: Gene object
        :param drug: Drug object
        :return: drug_gene_interaction_object or None
        """
        try:
            drug_gene_interaction_object = models.GeneDrugInteraction.objects.create(
                drug_id=drug,
                gene_id=gene,
                pubmed_id=pubmed_id,
                interaction_dataset_id=interaction_dataset_object
            )

            logging.debug('Created drug protein relation for drug {} and gene {} with id {}\n'.format(
                drug.id,
                gene.id,
                drug_gene_interaction_object.id
            ))
            return drug_gene_interaction_object

        except IntegrityError as e:
            # cancer type probably exists
            logging.debug(e)
            logging.debug('Could not create drug protein relation for drug {} and gene {}\n'.format(drug.id, gene.id))
            return

    def __check_and_add_cancer_dataset(self, dataset_name: str, link: str, version: str):
        """
        Interlude method
        Checks if cancer dataset entry exists in db and if so, returns it
        Otherwise creates new entry and returns it
        :param dataset_name:
        :return: cancer_dataset_object
        """

        # collect dataset id's in case they are missing
        cancer_dataset_object = self.__get_cancer_dataset(dataset_name)
        if cancer_dataset_object is None:
            # cancer dataset is not in db, add it
            cancer_dataset_object = self.__add_cancer_dataset(dataset_name, link, version)
            # if we dont have the cancer_dataset_object here, it is not in table but could also not be added
            assert isinstance(cancer_dataset_object, models.CancerDataset)

        return cancer_dataset_object

    def __check_and_add_cancer_types(self, unique_cancer_types: list, cancer_dataset_object: models.CancerDataset):
        """
        Interlude method
        Checks if cancer types are already in the db and if not,
        adds cancer types to the database by calling "__add_cancer_type".
        In any way, returns a lookup dict with {cancer_type: cancer_type_id}
        Returns
        :param unique_cancer_types: list of strings
        :return: lookup dict {cancer_type: cancer_type_object}
        """
        #
        # verify that cancer types are unique
        if len(unique_cancer_types) > len(set(unique_cancer_types)):
            raise ValueError('ERROR: unique_cancer_types contains duplicates')

        # collect cancer type id's
        cancer_types_objects = {}
        for cancer_type in unique_cancer_types:

            cancer_type_object = self.__get_cancer_type(cancer_type)
            if not cancer_type_object:
                # cancer type is not in db, add it
                cancer_type_object = self.__add_cancer_type(cancer_type)
            # if cancer_type_id is None, then something went wrong.
            # adding failed it is also not already in db
            assert isinstance(cancer_type_object, models.CancerType)
            cancer_types_objects[cancer_type] = cancer_type_object

            self.__add_cancer_type_dataset_relation(
                cancer_type=cancer_type_object,
                cancer_dataset=cancer_dataset_object
            )

        return cancer_types_objects

    def __check_and_add_cancer_genes(self,
                                     unique_cancer_genes: list,
                                     unique_entrez_ids: list):
        """
        Interlude method
        Checks if cancer genes are already in the db and if not,
        adds cancer genes to the database by calling "__add_gene".
        In any way, returns a lookup dict with {entrez_id: cancer_gene_id}
        Returns
        :param unique_cancer_genes: list of strings
        :param unique_entrez_ids: list of integers (pass list of None if no information)
        :return: lookup dict {entrez_id: cancer_gene_object}
        """

        # verify that cancer types are unique
        if len(unique_cancer_genes) > len(set(unique_cancer_genes)):
            raise ValueError('ERROR: unique_cancer_genes contains duplicates')

        # collect cancer entrez id's
        cancer_gene_objects = {}
        for index, cancer_gene in enumerate(unique_cancer_genes):
            entrez_id = unique_entrez_ids[index]

            cancer_gene_object = self.get_gene(entrez_id=entrez_id)

            if cancer_gene_object is None:
                # gene not yet in db, add it
                cancer_gene_object = self.__add_gene(gene=cancer_gene,
                                                     entrez_id=entrez_id)

            assert isinstance(cancer_gene_object, models.Gene)
            cancer_gene_objects[entrez_id] = cancer_gene_object

        return cancer_gene_objects

    def __check_and_add_drugs(self, drugs: list):
        """
        Interlude method
        gets a list of dictionaries that represent single drugs. each dictionary contains the all the drug attributes
        iterates over this list and tries to add the drugs to the database
        :param drugs:
        :return:
        """

        drug_objects = []
        for drug in drugs:

            drug_status_object, _ = models.DrugStatus.objects.get_or_create(
                name=drug['status']
            )
            if drug_status_object is None:
                return None
            try:
                pubchem_substance_id = int(drug['pubchem_substance_id'])
            except Exception:
                pubchem_substance_id = None

            drug_object, _ = models.Drug.objects.get_or_create(
                name=drug['name'],
                db_id=drug['drug_id'],
                status=drug_status_object,
                in_trial=drug['in_trial'],
                in_literature=drug['in_literature'],
                links=drug['links'],
                is_nutraceutical=drug['nutraceutical'],
                is_atc_antineoplastic_and_immunomodulating_agent=drug[
                    'is_atc_antineoplastic_and_immunomodulating_agent'
                ],
                pharmaGKD_id=drug['pharmaGKD_id'],
                pubchem_substance_id=pubchem_substance_id
            )

            drug_objects.append(drug_object)

        return drug_objects

    def add_gene_dataset(self, dataset_name, link, df, version):
        """
        Function to populate the database
        Add cancer genes given the dataset_name, a list of cancer genes and a list of cancer types

        :param df: a pandas dataframe with columns: name, cancer_type, entrez_id, pubmed_id
        :param link: reference to source
        :param dataset_name: string

        :return:
        """

        # create dataset name entry in case dataset name is new
        # collect dataset object
        cancer_dataset_object = self.__check_and_add_cancer_dataset(dataset_name, link, version)

        # replace all '_'  with ' '
        df['cancer_type'] = df['cancer_type'].map(lambda x: x.replace('_', ' '))
        unique_cancer_types = df['cancer_type'].unique()

        # create cancer types in db in case they are missing
        # collect cancer type objects
        cancer_types_objects = self.__check_and_add_cancer_types(unique_cancer_types, cancer_dataset_object)

        # add cancer genes
        # get unique entries
        unique_cancer_genes = []
        unique_entrez_ids = []
        # init entrez id as none since it is not mandatory, change if values are given
        for gene_name, entrez_id in zip(df['name'].to_list(), df['entrez_id'].to_list()):
            if gene_name in unique_cancer_genes:
                continue
            unique_cancer_genes.append(gene_name)

            unique_entrez_ids.append(entrez_id)

        unique_cancer_genes = df['name'].unique()
        cancer_genes_objects = self.__check_and_add_cancer_genes(unique_cancer_genes,
                                                                 unique_entrez_ids=unique_entrez_ids)

        # entries in cancer_genes and cancer_types belong together.
        # We have the lookup table from adding the cancer_types
        if 'pubmed_id' not in df.columns:
            df['pubmed_id'] = None
        for entrez_id, cancer_type, pubmed_id in list(zip(
                df['entrez_id'].to_list(), df['cancer_type'].to_list(), df['pubmed_id'].to_list())):
            cancer_type_object = cancer_types_objects[cancer_type]
            cancer_gene_object = cancer_genes_objects[entrez_id]

            # add gene entitys
            self.__add_cancer_gene_entity_instance(cancer_dataset_object, cancer_type_object, cancer_gene_object,
                                                   pubmed_id)

    def add_complete_drug_dataset(self, dataset_name, df):
        """
        Adds a complete drug dataset to the database table "drugs"
        Gets lists for each of the attributes where each index represents one drug

        :param df
        :return: drug_objects list(Drug)
        """

        # validate data
        # check for unique drug names
        if len(df['drug_id']) != len(df['drug_id'].unique()):
            raise ValueError('ERROR: Drug drug_id contains duplicates')

        # add drug dataset object
        drug_dataset_object, _ = models.DrugDataset.objects.get_or_create(name=dataset_name)

        # changing the format to a list of dictionaries
        drugs = df.to_dict(orient='records')

        drug_objects = self.__check_and_add_drugs(drugs)

        for drug_object in drug_objects:
            # establish connection from drug object to dataset object
            models.DrugEntity.objects.get_or_create(
                drug_id=drug_object,
                drug_dataset_id=drug_dataset_object
            )

        return

    def add_gene_gene_relations(self, df, interaction_dataset_name, link, version):
        """
        :param df: pandas Dataframe with columns: 'entrez_a','entrez_b','pubmed_id',
        :param interaction_dataset_name: str
        :param link: str, url to source
        confidence_value','detection_method_psi_mi','detection_method_name','type_psi_mi', 'type_name'
        :return:
        """
        # create entry for model InteractionGeneGeneDataset
        interaction_dataset_object, created = models.InteractionGeneGeneDataset.objects.get_or_create(
            name=interaction_dataset_name,
            link=link,
            version=version,
        )

        n = 0

        for index, row in df.iterrows():

            # TODO this is very hacky, change
            if 'from_protein_ac' in df.columns:
                gene_a = self.get_gene(uniprot_ac=row['from_protein_ac'])
                gene_b = self.get_gene(uniprot_ac=row['to_protein_ac'])

                if gene_a is None or gene_b is None or (gene_a == gene_b):
                    continue

                # check if a relation already exists
                gene_gene_interaction_object = self.__get_gene_gene_relation(gene_a, gene_b,
                                                                             interaction_dataset_object)

                if gene_gene_interaction_object is None:
                    # also check for the inverse relation
                    gene_gene_interaction_object = self.__get_gene_gene_relation(gene_b, gene_a,
                                                                                 interaction_dataset_object)

                if gene_gene_interaction_object is None:
                    gene_gene_interaction_object = self.__add_gene_gene_relation(
                        gene_a=gene_a,
                        gene_b=gene_b,
                        interaction_dataset_object=interaction_dataset_object
                    )
                    # assure that creation succeeded
                    assert isinstance(gene_gene_interaction_object, models.GeneGeneInteraction)
                    logging.debug('Created gene gene relation {}\n'.format(gene_gene_interaction_object.id))
                    n += 1
            elif 'gene_name_a' in df.columns:
                # this is the worst case since 'official' gene names are not really unique
                gene_a = models.Gene.objects.filter(name=row['gene_name_a'])
                gene_b = models.Gene.objects.filter(name=row['gene_name_b'])

                if len(gene_a) != 1 or len(gene_b) != 1 or (gene_a == gene_b):
                    continue
                gene_a = gene_a[0]
                gene_b = gene_b[0]

                # check if a relation already exists
                gene_gene_interaction_object = self.__get_gene_gene_relation(gene_a, gene_b,
                                                                             interaction_dataset_object)

                if gene_gene_interaction_object is None:
                    # also check for the inverse relation
                    gene_gene_interaction_object = self.__get_gene_gene_relation(gene_b, gene_a,
                                                                                 interaction_dataset_object)

                if gene_gene_interaction_object is None:
                    gene_gene_interaction_object = self.__add_gene_gene_relation(
                        gene_a=gene_a,
                        gene_b=gene_b,
                        interaction_dataset_object=interaction_dataset_object
                    )
                    # assure that creation succeeded
                    assert isinstance(gene_gene_interaction_object, models.GeneGeneInteraction)
                    logging.debug('Created gene gene relation {}\n'.format(gene_gene_interaction_object.id))
                    n += 1

            else:

                # skip self-interactions
                if row['entrez_a'] == row['entrez_b']:
                    continue

                gene_a = self.get_gene(entrez_id=row['entrez_a'])
                gene_b = self.get_gene(entrez_id=row['entrez_b'])

                # not all the genes must be in the database since we don't use the biogrid dataset to populate the genes
                if gene_a is None or gene_b is None:
                    continue

                # make sure we have found the cancer objects
                # if this fails, the protein is missing in the database
                assert isinstance(gene_a, models.Gene)
                assert isinstance(gene_b, models.Gene)

                # check if a relation already exists
                gene_gene_interaction_object = self.__get_gene_gene_relation(gene_a, gene_b,
                                                                             interaction_dataset_object)

                if gene_gene_interaction_object is None:
                    # also check for the inverse relation
                    gene_gene_interaction_object = self.__get_gene_gene_relation(gene_b, gene_a,
                                                                                 interaction_dataset_object)

                if gene_gene_interaction_object and 'pubmed_id' in row:
                    # if relation exists, add new values
                    gene_gene_interaction_object.pubmed_id += ',' + row['pubmed_id']
                    gene_gene_interaction_object.confidence_value += ',' + row['confidence_value']
                    gene_gene_interaction_object.detection_method_psi_mi += ',' + row['detection_method_psi_mi']
                    gene_gene_interaction_object.detection_method_name += ',' + row['detection_method_name']
                    gene_gene_interaction_object.type_psi_mi += ',' + row['type_psi_mi']
                    gene_gene_interaction_object.type_name += ',' + row['type_name']

                    # save the added data
                    gene_gene_interaction_object.save()

                    logging.debug('Data added to protein protein relation {}\n'.format(gene_gene_interaction_object.id))

                elif gene_gene_interaction_object is None:
                    # if there is no relation yet, create it
                    gene_gene_interaction_object = self.__add_gene_gene_relation(
                        gene_a=gene_a,
                        gene_b=gene_b,
                        interaction_dataset_object=interaction_dataset_object,
                        data=row
                    )

                    # assure that creation succeeded
                    assert isinstance(gene_gene_interaction_object, models.GeneGeneInteraction)
                    logging.debug('Created protein protein relation {}\n'.format(gene_gene_interaction_object.id))
                    n += 1
        interaction_dataset_object.n_interactions = n
        interaction_dataset_object.save()

    def add_drug_gene_relations_v2(self, df, interaction_dataset_name, link, version):
        """
        :param df: pandas Dataframe with column names 'drug_id', 'entrez_id'
        :return:
        """

        # create entry for model InteractionGeneGeneDataset
        interaction_dataset_object, created = models.InteractionGeneDrugDataset.objects.get_or_create(
            name=interaction_dataset_name,
            link=link,
            version=version
        )

        for index, row in df.iterrows():
            drug_object = self.__get_drug(db_id=row['drug_id'])
            gene_object = self.get_gene(entrez_id=row['entrez_id'])

            if gene_object is None or drug_object is None:
                continue

            # make sure we have found the objects
            # if this fails, the drug/protein is missing in the database
            assert isinstance(drug_object, models.Drug)
            assert isinstance(gene_object, models.Gene)

            drug_gene_interaction_object = self.__get_drug_gene_interaction(
                drug=drug_object,
                gene=gene_object,
                dataset=interaction_dataset_object
            )

            if drug_gene_interaction_object is None:
                # drug_gene_interaction_object does not exists yet
                drug_gene_interaction_object = self.__add_drug_gene_relation(
                    drug=drug_object,
                    gene=gene_object,
                    interaction_dataset_object=interaction_dataset_object
                )

                # assure that creation succeeded
                assert isinstance(drug_gene_interaction_object, models.GeneDrugInteraction)

    def add_drug_gene_relations(self, df, interaction_dataset_name, link, version, common=False):
        """
        :param df: pandas Dataframe with column names 'drug_id', 'entrez_id'
        :return:
        """

        # create entry for model InteractionGeneGeneDataset
        interaction_dataset_object, created = models.InteractionGeneDrugDataset.objects.get_or_create(
            name=interaction_dataset_name,
            link=link,
            version=version
        )

        for index, row in df.iterrows():

            if common:
                # TODO this is again suuuper hacky, change

                drug_object = self.__get_drug(db_id=row['drug_id'])
                gene_object = self.get_gene(uniprot_ac=row['protein_ac'])

                if gene_object is None or drug_object is None:
                    continue

                drug_gene_interaction_object = self.__get_drug_gene_interaction(
                    drug=drug_object,
                    gene=gene_object,
                    dataset=interaction_dataset_object
                )

                if drug_gene_interaction_object is None:
                    # drug_gene_interaction_object does not exists yet
                    drug_gene_interaction_object = self.__add_drug_gene_relation(
                        drug=drug_object,
                        gene=gene_object,
                        interaction_dataset_object=interaction_dataset_object
                    )

                    # assure that creation succeeded
                    assert isinstance(drug_gene_interaction_object, models.GeneDrugInteraction)

            else:

                drug_object = self.__get_drug(db_id=row['drug_id'])
                gene_object = self.get_gene(entrez_id=row['entrez_id'])

                if gene_object is None or drug_object is None:
                    continue

                # make sure we have found the objects
                # if this fails, the drug/protein is missing in the database
                assert isinstance(drug_object, models.Drug)
                assert isinstance(gene_object, models.Gene)

                drug_gene_interaction_object = self.__get_drug_gene_interaction(
                    drug=drug_object,
                    gene=gene_object,
                    dataset=interaction_dataset_object
                )

                if drug_gene_interaction_object is None:
                    # drug_gene_interaction_object does not exists yet
                    drug_gene_interaction_object = self.__add_drug_gene_relation(
                        drug=drug_object,
                        gene=gene_object,
                        interaction_dataset_object=interaction_dataset_object
                    )

                    # assure that creation succeeded
                    assert isinstance(drug_gene_interaction_object, models.GeneDrugInteraction)

    def add_common_gene_dataset(self, df):
        """
        Adds all information of df to Gene model

        :param df: pandas dataframe with columns gene_name, protein_ac, protein_name, entrez_id
        :return:
        """

        for index, row in df.iterrows():
            # add gene to db
            # here it is possible to add genes without entrez id (''). 
            if pd.isna(row['entrez_id']):
                row['entrez_id'] = None
            models.Gene.objects.get_or_create(
                name=row['gene_name'],
                entrez_id=row['entrez_id'],
                protein_name=row['protein_name'],
                uniprot_ac=row['protein_ac']
            )

    def add_comorbidity_disease_data(self, df):
        for index, row in df.iterrows():

            try:
                models.Disease.objects.create(
                    icd_10=row['ICD_10'],
                    mondo_id=row['Mondo'],
                    name=row['name']
                )
            except Exception as e:
                logging.debug(f'Could not create disease {row["name"]}. {e}')

    def add_comorbidity_disease_interactions(self, df):
        # 'disease1_ICD10', 'disease2_ICD10', 'phi_cor', 'relative_risk',
        #        'disease1_name', 'disease2_name', 'relative_risk_norm'
        for index, row in df.iterrows():
            disease_a = self.__get_disease(icd_10=row['disease1_ICD10'])
            disease_b = self.__get_disease(icd_10=row['disease2_ICD10'])

            # skip if one of the comorbidities is not in db
            if disease_a is None or disease_b is None:
                continue

            try:
                models.DiseaseInteractions.objects.create(
                    disease_a=disease_a,
                    disease_b=disease_b,
                    phi_cor=row['phi_cor'],
                    relative_risk=row['relative_risk'],
                    relative_risk_norm=row['relative_risk_norm']
                )
            except Exception as e:
                logging.debug(e)

    def add_comorbidity_gene_interactions(self, df):
        """
        dataframe has columns: 'EntrezID', 'Mondo', 'database_assertedBy', 'score_DisGeNet',
       'mondo_disease_name'
        :param df: pandas Dataframe
        :return:
        """
        for index, row in df.iterrows():
            # fetch gene
            gene_object = self.get_gene(entrez_id=row['Gene'])

            # fetch disease
            disease_object = self.__get_disease(mondo_id=row['Disease'])

            if disease_object is None:
                # create disease
                disease_object = models.Disease.objects.create(
                    mondo_id=row['Disease'],
                    name=row['Name']
                )

            # skip if gene or disease is not in db
            if gene_object is None or disease_object is None:
                continue
            try:
                models.DiseaseGeneInteractions.objects.create(
                    disease=disease_object,
                    gene=gene_object,
                    database_assertedBy=row['database_assertedBy'],
                    score_DisGeNet=row['score_DisGeNet']
                )
            except Exception as e:
                logging.debug(e)

    def add_shortest_distances(
            self,
            df,
            cancer_dataset_name,
            gene_interaction_dataset_name,
            drug_gene_interaction_dataset_name
    ):
        """

        :param df: columns source, target, dist, cancer_type
        :param cancer_dataset_name:
        :param gene_interaction_dataset_name:
        :param drug_gene_interaction_dataset_name:
        :return:
        """
        cancer_dataset = self.__get_cancer_dataset(cancer_dataset_name)
        gene_interaction_dataset = models.InteractionGeneGeneDataset.objects\
            .filter(name=gene_interaction_dataset_name)[0]
        drug_gene_interaction_dataset = models.InteractionGeneDrugDataset.objects\
            .filter(name=drug_gene_interaction_dataset_name)[0]

        for index, row in df.iterrows():
            # cancer_type might be none
            if pd.isnull(row['cancer_type']):
                continue
            else:
                cancer_type = models.CancerType.objects.get(id=row['cancer_type'])

            if row['source'][0] == 'g':
                # gene

                # graphIds are of format 'g<backendId>'
                gene_a = models.Gene.objects.get(id=row['source'][1:])
                gene_b = models.Gene.objects.get(id=row['target'][1:])

                models.ShortestDistanceGeneToCancerGene.objects.create(
                    cancer_dataset=cancer_dataset,
                    gene_interaction_dataset=gene_interaction_dataset,
                    cancer_type=cancer_type,
                    gene_a=gene_a,
                    gene_b=gene_b,
                    distance=row['dist']
                )

            else:
                # drug, 'd'
                # graphIds are of format 'g<backendId>'
                drug = models.Drug.objects.get(id=row['source'][1:])
                gene = models.Gene.objects.get(id=row['target'][1:])

                models.ShortestDistanceDrugToCancerGene.objects.create(
                    cancer_dataset=cancer_dataset,
                    gene_drug_interaction_dataset=drug_gene_interaction_dataset,
                    cancer_type=cancer_type,
                    drug=drug,
                    gene=gene,
                    distance=row['dist']
                )

    def add_expression(self, df_expression_norm, df_expression_tpm):
        def _add_expression_norm(df):
            for col in df:
                cancer_type, created = models.ExpressionCancerType.objects.get_or_create(
                    name=col
                )
                for gene_name, expr in df[col].iteritems():
                    for gene_model in models.Gene.objects.filter(name=gene_name).all():
                        try:
                            models.GeneExpressionLevel.objects.create(
                                gene=gene_model,
                                expression_cancer_type=cancer_type,
                                expression_level=expr
                            )
                        except IntegrityError:
                            pass

        def _add_expression_tpm(df):
            for col in df:
                cancer_type, created = models.ExpressionCancerType.objects.get_or_create(
                    name=col
                )
                for gene_name, expr in df[col].iteritems():
                    for gene_model in models.Gene.objects.filter(name=gene_name).all():
                        try:
                            ge = models.GeneExpressionLevel.objects.get(
                                gene=gene_model,
                                expression_cancer_type=cancer_type,
                            )
                            ge.expression_level_tpm = expr
                            ge.save()
                        except IntegrityError:
                            pass
        _add_expression_norm(df_expression_norm)
        _add_expression_tpm(df_expression_tpm)

    def add_tissues_and_expression(self, df):
        tissues_models = dict()
        for tissue_name in df.columns.values[2:]:
            try:
                tissue_model = models.Tissue.objects.get(name=tissue_name)
            except models.Tissue.DoesNotExist:
                tissue_model = models.Tissue.objects.create(name=tissue_name)
            tissues_models[tissue_name] = tissue_model

        for _, row in df.iterrows():
            gene_name = row['UniprotId']
            for gene_model in models.Gene.objects.filter(name=gene_name).all():
                for tissue_name, tissue_model in tissues_models.items():
                    try:
                        models.ExpressionLevel.objects.create(
                            gene=gene_model,
                            tissue=tissue_model,
                            expression_level=row[tissue_name]
                        )
                    except IntegrityError:
                        pass

    def add_mutations_counts(self, df):
        for index, row in df.iterrows():
            mutation_cancer_type_object, created = models.MutationCancerType.objects.get_or_create(
                name=row['name'],
                abbreviation=row['cancer.type']
            )
            gene = self.get_gene(entrez_id=row['entrez_id'])
            # check if gene was found, we might not have all genes in our db
            if gene is None:
                continue

            count, count_created = models.MutationCounts.objects.get_or_create(
                cancer_type=mutation_cancer_type_object,
                gene=gene,
                mutation_counts=row['mutation_count'],
                mutation_score=row['mutation_score']
            )

    def add_ctrpv2_data(self, df):
        for _, row in df.iterrows():
            drug_objects = models.Drug.objects.filter(name__iexact=row['cpd_name'])
            if drug_objects:
                print(row['cpd_name'])
                drug_object = drug_objects[0]
                drug_object.ctrpv2_id = int(row['master_cpd_id'])
                drug_object.save()

    def add_cancernet_table(self, df, region, targeted):
        for index, row in df.iterrows():
            # find drug
            drug_object = None
            try:
                drug_object = models.Drug.objects.get(name=row['Drug name'])
            except models.Drug.DoesNotExist:
                continue
            # if not drug_object:
            #     try:
            #         drug_object = models.Drug.objects.get(name=row['Drug name'].split(' ')[0])
            #     except models.Drug.DoesNotExist:
            #         pass
            models.Cancernet.objects.create(
                drug=drug_object,
                target=row['Target '],
                cancer_type_long=row['Cancer type (from website)'],
                cancer_type=row['Cancer type (short version)'],
                link=row['website'],
                access_date=row['Access date'],
                approved_region=region,
                notes=row['Notes'],
                combination_formation=row['Combination formation'],
                targeted=targeted
            )

    def add_drug_gene_actions(self, df, dataset):
        dataset_object = models.InteractionGeneDrugDataset.objects.get(name=dataset)
        for _, row in df.iterrows():
            try:
                gene = models.Gene.objects.get(uniprot_ac=row['uniprot_id'])
                drug = models.Drug.objects.get(db_id=row['drugbank_id'])
            except (models.Gene.DoesNotExist, models.Drug.DoesNotExist):
                continue
            try:
                interaction_object = models.GeneDrugInteraction.objects.get(
                    gene_id=gene,
                    drug_id=drug,
                    interaction_dataset_id=dataset_object
                )
                interaction_object.action = row['actions']
                interaction_object.save()
            except models.GeneDrugInteraction.DoesNotExist:
                continue

    @staticmethod
    def normalize_tissue_expression():
        for t in models.Tissue.objects.all():
            exps = models.ExpressionLevel.objects.filter(tissue=t)
            exps_lvls = [np.log2(x.expression_level+1) if x.expression_level else 0 for x in exps]
            exps_lvls_norm = np.array(exps_lvls) / max(exps_lvls)
            for i, e in enumerate(exps_lvls_norm):
                exps[i].expression_level_norm = e
                exps[i].save()

    @staticmethod
    def internal_drugs_all():
        # fetch drugs and drug interactions
        drug_interaction_objects = models.GeneDrugInteraction.objects.all()
        drug_interactions = serializers.GeneDrugInteractionSerializer(many=True) \
            .to_representation(drug_interaction_objects)
        """
        # add interaction type
        for drug_interaction in drug_interactions:
            # drug-node type
            drug_interaction['type'] = 'd-n'
        """
        drug_objects = models.Drug.objects.all()
        drugs = serializers.DrugSerializer(many=True).to_representation(drug_objects)

        return {
            'drugs': drugs,
            'drugEdges': drug_interactions
        }

    @staticmethod
    def internal_drugs(drug_interaction_dataset_name):
        # fetch drugs and drug interactions
        drug_interaction_objects = models.GeneDrugInteraction.objects.filter(
            interaction_dataset_id__name=drug_interaction_dataset_name
        )
        drug_interactions = serializers.GeneDrugInteractionSerializer(many=True) \
            .to_representation(drug_interaction_objects)
        """
        # add interaction type
        for drug_interaction in drug_interactions:
            # drug-node type
            drug_interaction['type'] = 'd-n'
        """
        drug_objects = models.Drug.objects.all()
        drugs = serializers.DrugSerializer(many=True).to_representation(drug_objects)

        return {
            'drugs': drugs,
            'drugEdges': drug_interactions
        }

    @staticmethod
    def internal_genes(gene_dataset_name):
        # fetch all genes for cancer gene dataset
        cancer_node_objects = models.CancerGeneEntity.objects.filter(cancer_dataset_id__name=gene_dataset_name)
        cancer_nodes = serializers.CancerGeneEntitySerializer(many=True).to_representation(cancer_node_objects)
        cancer_node_ids = [gene['backendId'] for gene in cancer_nodes]

        # get all node objects but exclude the cancer nodes
        node_objects = models.Gene.objects.exclude(id__in=cancer_node_ids)
        nodes = serializers.GeneSerializer(many=True).to_representation(node_objects)

        response = {
            'nodes': nodes,
            'cancerNodes': cancer_nodes,
        }
        return response

    @staticmethod
    def internal_gene_interactions_all():
        # get all interactions
        node_node_interaction_objects = models.GeneGeneInteraction.objects.order_by('gene_b')
        node_node_interactions = serializers.GeneGeneInteractionSerializer(many=True) \
            .to_representation(node_node_interaction_objects)
        # merge dataset names
        df = pd.DataFrame.from_records(node_node_interactions)
        # df['dataset_name_internal'] = df[
        #     ['interactor_a_graphId', 'interactor_b_graphId', 'dataset_name_internal']
        #     ].groupby(
        #         ['interactor_a_graphId', 'interactor_b_graphId']
        #         )['dataset_name_internal'].transform(lambda x: ','.join(x))
        df['dataset_name'] = df[
            ['interactor_a_graphId', 'interactor_b_graphId', 'dataset_name']
            ].groupby(
                ['interactor_a_graphId', 'interactor_b_graphId']
                )['dataset_name'].transform(lambda x: ','.join(x))
        df = df.drop_duplicates(subset=['interactor_a_graphId', 'interactor_b_graphId'])
        return df.to_dict('records')

    @staticmethod
    def internal_gene_interactions(gene_interaction_dataset_name):
        # TODO put this function into extra file
        # get all interactions
        node_node_interaction_objects = models.GeneGeneInteraction.objects.filter(
            interaction_dataset_id__name=gene_interaction_dataset_name
        )
        node_node_interactions = serializers.GeneGeneInteractionSerializer(many=True) \
            .to_representation(node_node_interaction_objects)

        return {'edges': node_node_interactions}

    @staticmethod
    def internal_gene_scores(gene_backendId):
        gene_object = models.Gene.objects.get(id=int(gene_backendId))

        tissue_objects = models.Tissue.objects.all()
        tissue_expression_scores = {t.name: None for t in tissue_objects}
        for t in tissue_objects:
            res = models.ExpressionLevel.objects.filter(
                tissue=t,
                gene=gene_object
            )
            if res:
                tissue_expression_scores[t.name] = res[0].expression_level_norm

        # get expression scores
        expression_cancer_type_objects = models.ExpressionCancerType.objects.all()
        expression_cancer_scores = {t.name: None for t in expression_cancer_type_objects}
        for t in expression_cancer_type_objects:
            res = models.GeneExpressionLevel.objects.filter(
                expression_cancer_type=t,
                gene=gene_object
            )
            if res:
                expression_cancer_scores[t.name] = res[0].expression_level

        # get mutation scores
        mutation_type_objects = models.MutationCancerType.objects.all()
        mutation_scores = {mt.name: None for mt in mutation_type_objects}
        for mt in mutation_type_objects:
            res = models.MutationCounts.objects.filter(
                cancer_type=mt,
                gene=gene_object
            )
            if res:
                mutation_scores[mt.name] = res[0].mutation_score

        return {
            'mutation_scores': mutation_scores,
            'cancer_expression_scores': expression_cancer_scores,
            'expression_scores': tissue_expression_scores
            }
