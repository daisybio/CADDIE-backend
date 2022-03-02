# Serializers define the API representation.
import json

from rest_framework import serializers

import caddie.models as models


class CancerDatasetSerializer(serializers.ModelSerializer):
    backendId = serializers.SerializerMethodField()
    name = serializers.SerializerMethodField()
    link = serializers.SerializerMethodField()
    version = serializers.SerializerMethodField()

    def get_backendId(self, obj):
        return obj.id

    def get_name(self, obj):
        return obj.name

    def get_link(self, obj):
        return obj.link

    def get_version(self, obj):
        return obj.version

    class Meta:
        model = models.CancerDataset
        fields = ['backendId', 'name', 'link', 'version']


class InteractionDatasetSerializer(serializers.ModelSerializer):
    backendId = serializers.SerializerMethodField()
    name = serializers.SerializerMethodField()
    link = serializers.SerializerMethodField()
    version = serializers.SerializerMethodField()

    def get_backendId(self, obj):
        return obj.id

    def get_name(self, obj):
        return obj.name

    def get_link(self, obj):
        return obj.link

    def get_version(self, obj):
        return obj.version

    class Meta:
        model = models.InteractionGeneGeneDataset
        fields = ['backendId', 'name', 'link', 'version']


class GeneSerializer(serializers.ModelSerializer):
    backendId = serializers.SerializerMethodField()
    name = serializers.SerializerMethodField()
    entrez_id = serializers.SerializerMethodField()
    gene_alias = serializers.SerializerMethodField()
    protein_name = serializers.SerializerMethodField()
    uniprot_ac = serializers.SerializerMethodField()
    graphId = serializers.SerializerMethodField()

    def get_backendId(self, obj):
        return obj.id

    def get_graphId(self, obj):
        return f'g{obj.id}'

    def get_name(self, obj):
        return obj.name

    def get_entrez_id(self, obj):
        return obj.entrez_id

    def get_gene_alias(self, obj):
        return obj.alias

    def get_protein_name(self, obj):
        return obj.protein_name

    def get_uniprot_ac(self, obj):
        return obj.uniprot_ac

    class Meta:
        model = models.Gene
        fields = ['backendId', 'graphId', 'name', 'entrez_id', 'gene_alias', 'protein_name', 'uniprot_ac']


class DiseaseSerializer(serializers.ModelSerializer):
    backendId = serializers.SerializerMethodField()
    name = serializers.SerializerMethodField()
    icd_10 = serializers.SerializerMethodField()
    mondo_id = serializers.SerializerMethodField()

    def get_backendId(self, obj):
        return obj.id

    def get_name(self, obj):
        return obj.name

    def get_icd_10(self, obj):
        return obj.icd_10

    def get_mondo_id(self, obj):
        return obj.mondo_id

    class Meta:
        model = models.Disease
        fields = ['backendId', 'name', 'icd_10', 'mondo_id']


class DiseaseGeneInteractionsSerializer(serializers.ModelSerializer):
    disease_backendId = serializers.SerializerMethodField()
    disease_name = serializers.SerializerMethodField()
    disease_mondo_id = serializers.SerializerMethodField()
    disease_icd_10 = serializers.SerializerMethodField()
    gene_entrez_id = serializers.SerializerMethodField()
    gene_name = serializers.SerializerMethodField()
    gene_backendId = serializers.SerializerMethodField()

    def get_gene_backendId(self, obj):
        return obj.gene.id

    def get_gene_name(self, obj):
        return obj.gene.name

    def get_gene_entrez_id(self, obj):
        return obj.gene.entrez_id

    def get_disease_backendId(self, obj):
        return obj.disease.id

    def get_disease_name(self, obj):
        return obj.disease.name

    def get_disease_mondo_id(self, obj):
        return obj.disease.mondo_id

    def get_disease_icd_10(self, obj):
        return obj.disease.icd_10

    class Meta:
        model = models.DiseaseGeneInteractions
        fields = ['disease_backendId', 'disease_name', 'disease_mondo_id', 'disease_icd_10',
                  'gene_entrez_id', 'gene_name', 'gene_backendId']


class CancerGeneEntitySerializer(serializers.ModelSerializer):
    backendId = serializers.SerializerMethodField()
    graphId = serializers.SerializerMethodField()
    name = serializers.SerializerMethodField()
    entrez_id = serializers.SerializerMethodField()
    gene_alias = serializers.SerializerMethodField()
    type = serializers.SerializerMethodField()
    type_id = serializers.SerializerMethodField()
    cancer_occurrences = serializers.SerializerMethodField()
    protein_name = serializers.SerializerMethodField()
    uniprot_ac = serializers.SerializerMethodField()

    def get_backendId(self, obj):
        return obj.cancer_gene_id.id

    def get_graphId(self, obj):
        return f'g{obj.cancer_gene_id.id}'

    def get_name(self, obj):
        return obj.cancer_gene_id.name

    def get_entrez_id(self, obj):
        return obj.cancer_gene_id.entrez_id

    def get_gene_alias(self, obj):
        return obj.cancer_gene_id.alias

    def get_type(self, obj):
        return obj.cancer_type_id.name

    def get_type_id(self, obj):
        return obj.cancer_type_id.id

    def get_cancer_occurrences(self, obj):
        return obj.cancer_occurrences

    def get_uniprot_ac(self, obj):
        return obj.cancer_gene_id.uniprot_ac

    def get_protein_name(self, obj):
        return obj.cancer_gene_id.protein_name

    class Meta:
        model = models.CancerGeneEntity
        fields = [
            'backendId',
            'graphId',
            'name',
            'entrez_id',
            'gene_alias',
            'type',
            'type_id',
            'cancer_occurrences',
            'uniprot_ac',
            'protein_name'
        ]


class RelationDatasetTypeSerializer(serializers.ModelSerializer):
    backendId = serializers.SerializerMethodField()
    name = serializers.SerializerMethodField()

    def get_backendId(self, obj):
        return obj.cancer_type_id.id

    def get_name(self, obj):
        return obj.cancer_type_id.name

    class Meta:
        model = models.RelationDatasetType
        fields = ['backendId', 'name']


class CancernetSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.Cancernet
        fields = '__all__'


class CancerTypeSerializer(serializers.ModelSerializer):
    backendId = serializers.SerializerMethodField()
    name = serializers.SerializerMethodField()

    def get_backendId(self, obj):
        return obj.id

    def get_name(self, obj):
        return obj.name

    class Meta:
        model = models.CancerType
        fields = ['backendId', 'name']


class DrugSerializer(serializers.ModelSerializer):
    backendId = serializers.SerializerMethodField()
    db_id = serializers.SerializerMethodField()
    name = serializers.SerializerMethodField()
    status = serializers.SerializerMethodField()
    in_trial = serializers.SerializerMethodField()
    in_literature = serializers.SerializerMethodField()
    links = serializers.SerializerMethodField()
    graphId = serializers.SerializerMethodField()
    in_cancernet = serializers.SerializerMethodField()
    is_nutraceutical = serializers.SerializerMethodField()
    is_atc_antineoplastic_and_immunomodulating_agent = serializers.SerializerMethodField()

    def get_backendId(self, obj):
        return obj.id

    def get_graphId(self, obj):
        return f'd{obj.id}'

    def get_db_id(self, obj):
        return obj.db_id

    def get_name(self, obj):
        return obj.name

    def get_status(self, obj):
        return obj.status.name

    def get_in_trial(self, obj):
        return obj.in_trial

    def get_in_literature(self, obj):
        return obj.in_literature

    def get_links(self, obj):
        return obj.links

    def get_in_cancernet(self, obj):
        return len(models.Cancernet.objects.filter(drug=int(obj.id))) > 0

    def get_is_nutraceutical(self, obj):
        return obj.is_nutraceutical

    def get_is_atc_antineoplastic_and_immunomodulating_agent(self, obj):
        return obj.is_atc_antineoplastic_and_immunomodulating_agent

    class Meta:
        model = models.Drug
        fields = [
            'backendId',
            'graphId',
            'db_id',
            'name',
            'status',
            'in_trial',
            'in_literature',
            'links',
            'in_cancernet',
            'is_nutraceutical',
            'is_atc_antineoplastic_and_immunomodulating_agent'
        ]


class MutationCancerTypeSerializer(serializers.ModelSerializer):
    backendId = serializers.SerializerMethodField()

    def get_backendId(self, obj):
        return obj.id

    def get_abbreviation(self, obj):
        return obj.abbreviation

    class Meta:
        model = models.MutationCancerType
        fields = ['backendId', 'name', 'abbreviation']


class MutationCountsGeneSerializer(serializers.ModelSerializer):
    backendId = serializers.SerializerMethodField()
    name = serializers.SerializerMethodField()
    entrez_id = serializers.SerializerMethodField()
    protein_name = serializers.SerializerMethodField()
    uniprot_ac = serializers.SerializerMethodField()
    graphId = serializers.SerializerMethodField()
    mutation_counts = serializers.SerializerMethodField()
    mutation_score = serializers.SerializerMethodField()

    def get_backendId(self, obj):
        return obj.gene.id

    def get_name(self, obj):
        return obj.gene.name

    def get_entrez_id(self, obj):
        return obj.gene.entrez_id

    def get_protein_name(self, obj):
        return obj.gene.protein_name

    def get_uniprot_ac(self, obj):
        return obj.gene.uniprot_ac

    def get_graphId(self, obj):
        return f'g{obj.gene.id}'

    def get_mutation_counts(self, obj):
        return obj.mutation_counts

    def get_mutation_score(self, obj):
        return obj.mutation_score

    class Meta:
        model = models.MutationCancerType
        fields = [
            'backendId',
            'name',
            'entrez_id',
            'protein_name',
            'uniprot_ac',
            'graphId',
            'mutation_counts',
            'mutation_score'
        ]


class TissueSerializer(serializers.ModelSerializer):
    backendId = serializers.SerializerMethodField()

    def get_backendId(self, obj):
        return obj.id

    class Meta:
        model = models.Tissue
        fields = ['backendId', 'name']


class DrugTargetActionSerializer(serializers.ModelSerializer):
    backendId = serializers.SerializerMethodField()

    def get_backendId(self, obj):
        return obj.id

    class Meta:
        model = models.DrugTargetAction
        fields = ['backendId', 'name']


class ExpressionCancerTypeSerializer(serializers.ModelSerializer):
    backendId = serializers.SerializerMethodField()

    def get_backendId(self, obj):
        return obj.id

    class Meta:
        model = models.ExpressionCancerType
        fields = ['backendId', 'name']


class DrugStatusSerializer(serializers.ModelSerializer):
    backendId = serializers.SerializerMethodField()

    def get_backendId(self, obj):
        return obj.id

    class Meta:
        model = models.Tissue
        fields = ['backendId', 'name']


class DrugEntitySerializer(serializers.ModelSerializer):
    backendId = serializers.SerializerMethodField()
    db_id = serializers.SerializerMethodField()
    name = serializers.SerializerMethodField()
    status = serializers.SerializerMethodField()
    in_trial = serializers.SerializerMethodField()
    in_literature = serializers.SerializerMethodField()
    in_cancernet = serializers.SerializerMethodField()
    links = serializers.SerializerMethodField()

    def get_backendId(self, obj):
        return obj.drug_id.id

    def get_db_id(self, obj):
        return obj.drug_id.db_id

    def get_name(self, obj):
        return obj.drug_id.name

    def get_status(self, obj):
        return obj.drug_id.status.name

    def get_in_trial(self, obj):
        return obj.drug_id.in_trial

    def get_in_literature(self, obj):
        return obj.drug_id.in_literature

    def get_in_cancernet(self, obj):
        return len(models.Cancernet.objects.filter(drug=int(obj.drug_id))) > 0

    def get_links(self, obj):
        return obj.drug_id.links

    class Meta:
        model = models.DrugEntity
        fields = ['backendId', 'db_id', 'name', 'status', 'in_trial', 'in_literature', 'links', 'in_cancernet']


class ShortestDistanceDrugToCancerGeneSerializer(serializers.ModelSerializer):
    cancer_dataset_name = serializers.SerializerMethodField()
    cancer_dataset_backendId = serializers.SerializerMethodField()

    gene_drug_interaction_dataset_name = serializers.SerializerMethodField()
    gene_drug_interaction_dataset_backendId = serializers.SerializerMethodField()

    cancer_type_name = serializers.SerializerMethodField()
    cancer_type_backendId = serializers.SerializerMethodField()

    gene_entrezId = serializers.SerializerMethodField()
    gene_graphId = serializers.SerializerMethodField()
    gene_name = serializers.SerializerMethodField()
    gene_backendId = serializers.SerializerMethodField()

    drug_db_id = serializers.SerializerMethodField()
    drug_graphId = serializers.SerializerMethodField()
    drug_name = serializers.SerializerMethodField()
    drug_backendId = serializers.SerializerMethodField()

    distance = serializers.SerializerMethodField()

    def get_cancer_dataset_backendId(self, obj):
        return obj.cancer_dataset.name

    def get_cancer_dataset_name(self, obj):
        return obj.cancer_dataset.id

    def get_gene_drug_interaction_dataset_name(self, obj):
        return obj.gene_drug_interaction_dataset.name

    def get_gene_drug_interaction_dataset_backendId(self, obj):
        return obj.gene_drug_interaction_dataset.id

    def get_cancer_type_name(self, obj):
        return obj.cancer_type.name

    def get_cancer_type_backendId(self, obj):
        return obj.cancer_type.id

    def get_gene_graphId(self, obj):
        return f'g{obj.gene.id}'

    def get_gene_entrezId(self, obj):
        return obj.gene.entrez_id

    def get_gene_name(self, obj):
        return obj.gene.name

    def get_gene_backendId(self, obj):
        return obj.gene.id

    def get_drug_graphId(self, obj):
        return f'd{obj.drug.id}'

    def get_drug_db_id(self, obj):
        return obj.drug.db_id

    def get_drug_name(self, obj):
        return obj.drug.name

    def get_drug_backendId(self, obj):
        return obj.drug.id

    def get_distance(self, obj):
        return obj.distance

    class Meta:
        model = models.ShortestDistanceGeneToCancerGene
        fields = [
            'cancer_dataset_name',
            'cancer_dataset_backendId',

            'gene_drug_interaction_dataset_name',
            'gene_drug_interaction_dataset_backendId',

            'cancer_type_name',
            'cancer_type_backendId',

            'gene_entrezId',
            'gene_graphId',
            'gene_name',
            'gene_backendId',

            'drug_db_id',
            'drug_graphId',
            'drug_name',
            'drug_backendId',

            'distance'
        ]


class ShortestDistanceCancerGeneSerializer(serializers.ModelSerializer):
    cancer_dataset_name = serializers.SerializerMethodField()
    cancer_dataset_backendId = serializers.SerializerMethodField()

    gene_interaction_dataset_name = serializers.SerializerMethodField()
    gene_interaction_dataset_backendId = serializers.SerializerMethodField()

    cancer_type_name = serializers.SerializerMethodField()
    cancer_type_backendId = serializers.SerializerMethodField()

    gene_a_entrezId = serializers.SerializerMethodField()
    gene_a_graphId = serializers.SerializerMethodField()
    gene_a_name = serializers.SerializerMethodField()
    gene_a_backendId = serializers.SerializerMethodField()

    gene_b_entrezId = serializers.SerializerMethodField()
    gene_b_graphId = serializers.SerializerMethodField()
    gene_b_name = serializers.SerializerMethodField()
    gene_b_backendId = serializers.SerializerMethodField()

    distance = serializers.SerializerMethodField()

    def get_cancer_dataset_backendId(self, obj):
        return obj.cancer_dataset.name

    def get_cancer_dataset_name(self, obj):
        return obj.cancer_dataset.id

    def get_gene_interaction_dataset_name(self, obj):
        return obj.gene_interaction_dataset.name

    def get_gene_interaction_dataset_backendId(self, obj):
        return obj.gene_interaction_dataset.id

    def get_cancer_type_name(self, obj):
        return obj.cancer_type.name

    def get_cancer_type_backendId(self, obj):
        return obj.cancer_type.id

    def get_gene_a_graphId(self, obj):
        return f'g{obj.gene_a.id}'

    def get_gene_a_entrezId(self, obj):
        return obj.gene_a.entrez_id

    def get_gene_a_name(self, obj):
        return obj.gene_a.name

    def get_gene_a_backendId(self, obj):
        return obj.gene_a.id

    def get_gene_b_graphId(self, obj):
        return f'g{obj.gene_b.id}'

    def get_gene_b_entrezId(self, obj):
        return obj.gene_b.entrez_id

    def get_gene_b_name(self, obj):
        return obj.gene_b.name

    def get_gene_b_backendId(self, obj):
        return obj.gene_b.id

    def get_distance(self, obj):
        return obj.distance

    class Meta:
        model = models.ShortestDistanceGeneToCancerGene
        fields = [
            'cancer_dataset_name',
            'cancer_dataset_backendId',

            'gene_interaction_dataset_name',
            'gene_interaction_dataset_backendId',

            'cancer_type_name',
            'cancer_type_backendId',

            'gene_a_entrezId',
            'gene_a_graphId',
            'gene_a_name',
            'gene_a_backendId',

            'gene_b_entrezId',
            'gene_b_graphId',
            'gene_b_name',
            'gene_b_backendId',

            'distance'
        ]


class GeneGeneInteractionSerializer(serializers.ModelSerializer):
    interactor_a_name = serializers.SerializerMethodField()
    interactor_a_backendId = serializers.SerializerMethodField()
    interactor_a_entrez_id = serializers.SerializerMethodField()
    interactor_a_graphId = serializers.SerializerMethodField()
    interactor_a_uniprot_ac = serializers.SerializerMethodField()
    interactor_a_protein_name = serializers.SerializerMethodField()
    interactor_b_name = serializers.SerializerMethodField()
    interactor_b_backendId = serializers.SerializerMethodField()
    interactor_b_entrez_id = serializers.SerializerMethodField()
    interactor_b_graphId = serializers.SerializerMethodField()
    interactor_b_uniprot_ac = serializers.SerializerMethodField()
    interactor_b_protein_name = serializers.SerializerMethodField()
    backendId = serializers.SerializerMethodField()
    graphId = serializers.SerializerMethodField()
    dataset_name = serializers.SerializerMethodField()

    def get_backendId(self, obj):
        return obj.id

    def get_graphId(self, obj):
        return f'gg{obj.id}'

    def get_interactor_a_name(self, obj):
        return obj.gene_a.name

    def get_interactor_b_name(self, obj):
        return obj.gene_b.name

    def get_interactor_a_backendId(self, obj):
        return obj.gene_a.id

    def get_interactor_b_backendId(self, obj):
        return obj.gene_b.id

    def get_interactor_a_entrez_id(self, obj):
        return obj.gene_a.entrez_id

    def get_interactor_b_entrez_id(self, obj):
        return obj.gene_b.entrez_id

    def get_interactor_a_graphId(self, obj):
        return f'g{obj.gene_a.id}'

    def get_interactor_b_graphId(self, obj):
        return f'g{obj.gene_b.id}'

    def get_interactor_a_uniprot_ac(self, obj):
        return obj.gene_a.uniprot_ac

    def get_interactor_b_uniprot_ac(self, obj):
        return obj.gene_b.uniprot_ac

    def get_interactor_a_protein_name(self, obj):
        return obj.gene_a.protein_name

    def get_interactor_b_protein_name(self, obj):
        return obj.gene_b.protein_name

    def get_dataset_name(self, obj):
        return obj.interaction_dataset_id.name

    class Meta:
        model = models.GeneGeneInteraction
        fields = [
            'interactor_a_name',
            'interactor_a_backendId',
            'interactor_a_entrez_id',
            'interactor_a_graphId',
            'interactor_a_uniprot_ac',
            'interactor_a_protein_name',
            'interactor_b_name',
            'interactor_b_backendId',
            'interactor_b_entrez_id',
            'interactor_b_graphId',
            'interactor_b_uniprot_ac',
            'interactor_b_protein_name',
            'backendId',
            'graphId',
            'dataset_name'
            ]


class GeneDrugInteractionSerializer(serializers.ModelSerializer):
    gene_entrezId = serializers.SerializerMethodField()
    drug_dbId = serializers.SerializerMethodField()
    drug_backendId = serializers.SerializerMethodField()
    gene_backendId = serializers.SerializerMethodField()
    pubmed_id = serializers.SerializerMethodField()
    gene_graphId = serializers.SerializerMethodField()
    drug_graphId = serializers.SerializerMethodField()
    backendId = serializers.SerializerMethodField()
    graphId = serializers.SerializerMethodField()
    dataset_name = serializers.SerializerMethodField()

    def get_backendId(self, obj):
        return obj.id

    def get_graphId(self, obj):
        return f'dg{obj.id}'

    def get_gene_entrezId(self, obj):
        return ''  # obj.gene_id.entrez_id

    def get_drug_dbId(self, obj):
        return obj.drug_id.db_id

    def get_drug_backendId(self, obj):
        return obj.drug_id.id

    def get_drug_graphId(self, obj):
        return f'd{obj.drug_id.id}'

    def get_gene_backendId(self, obj):
        return obj.gene_id.id

    def get_gene_graphId(self, obj):
        return f'g{obj.gene_id.id}'

    def get_pubmed_id(self, obj):
        return ''   # obj.pubmed_id

    def get_dataset_name(self, obj):
        return obj.interaction_dataset_id.name

    class Meta:
        model = models.GeneDrugInteraction
        fields = [
            'gene_entrezId',
            'drug_dbId',
            'drug_backendId',
            'gene_backendId',
            'pubmed_id',
            'gene_graphId',
            'drug_graphId',
            'backendId',
            'graphId',
            'action',
            'dataset_name'
        ]


class RelationGeneDrugSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.GeneDrugInteraction
        fields = '__all__'


class TaskSerializer(serializers.ModelSerializer):
    parameters = serializers.SerializerMethodField()

    def get_parameters(self, obj):
        parameters = json.loads(obj.parameters)
        # add names for parameters
        if 'source_tasks' in parameters:
            return parameters
        cancer_type_names = [t.name for t in models.CancerType.objects.filter(id__in=parameters['cancer_types'])]
        parameters['cancer_type_names'] = cancer_type_names
        return parameters

    class Meta:
        model = models.Task
        fields = ['algorithm', 'target', 'parameters', 'job_id', 'worker_id', 'progress', 'status', 'created_at',
                  'started_at', 'finished_at', 'done', 'failed']


class TaskStatusSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.Task
        fields = ['algorithm', 'target', 'progress', 'status', 'created_at', 'started_at', 'finished_at', 'done',
                  'failed']
