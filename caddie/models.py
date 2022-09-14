from django.db import models
from django.contrib.contenttypes.fields import GenericForeignKey
from django.contrib.contenttypes.models import ContentType


class DrugDataset(models.Model):
    name = models.CharField(max_length=128, default='', unique=True)


class DrugStatus(models.Model):
    name = models.CharField(max_length=128, default='', unique=True)


class Drug(models.Model):
    db_id = models.CharField(max_length=7, default='', unique=True)
    name = models.CharField(max_length=1000, default='')
    status = models.ForeignKey('DrugStatus', on_delete=models.CASCADE)
    in_trial = models.BooleanField(default=False)
    in_literature = models.BooleanField(default=False)
    links = models.CharField(max_length=5000, default='')
    is_nutraceutical = models.BooleanField(default=False)
    is_atc_antineoplastic_and_immunomodulating_agent = models.BooleanField(default=False)
    pharmaGKD_id = models.CharField(max_length=20, default='', unique=False)
    pubchem_substance_id = models.PositiveIntegerField(default=None, unique=True, null=True, blank=True)
    ctrpv2_id = models.PositiveIntegerField(default=None, null=True)


class DrugEntity(models.Model):
    drug_id = models.ForeignKey('Drug', on_delete=models.CASCADE, related_name='drug_id')
    drug_dataset_id = models.ForeignKey('DrugDataset', on_delete=models.CASCADE, related_name='drug_dataset_id')


class Cancernet(models.Model):
    drug = models.ForeignKey('Drug', on_delete=models.CASCADE, related_name='drug')
    target = models.CharField(max_length=1000, default='')
    cancer_type_long = models.CharField(max_length=1000, default='')
    cancer_type = models.CharField(max_length=100, default='')
    link = models.CharField(max_length=500, default='')
    access_date = models.CharField(max_length=10, default='')
    approved_region = models.CharField(max_length=10, default='')  # NIH, UK
    notes = models.CharField(max_length=1000, default='')
    combination_formation = models.CharField(max_length=500, default='')
    targeted = models.BooleanField(default=False)


class CancerDataset(models.Model):
    name = models.CharField(max_length=128, default='', unique=True)
    link = models.CharField(max_length=128, default='', unique=True)
    version = models.CharField(max_length=64, default='')


class CancerType(models.Model):
    name = models.CharField(max_length=128, default='', unique=True)


class Gene(models.Model):
    name = models.CharField(max_length=128, default='', unique=False)
    entrez_id = models.PositiveIntegerField(default=None, unique=True, null=True, blank=True)
    alias = models.CharField(max_length=2048, default='')
    protein_name = models.CharField(max_length=128, default='', unique=False, null=True)
    uniprot_ac = models.CharField(max_length=128, default='', unique=True, null=True)


class Tissue(models.Model):
    """
    Tissue objects for expression per tissue information
    """
    name = models.CharField(max_length=128, default='')


class ExpressionCancerType(models.Model):
    """
    Cancer Type objects for expression per tissue information
    """
    name = models.CharField(max_length=128, default='')
    abbreviation = models.CharField(max_length=128, default='')


class ExpressionLevel(models.Model):
    """
    Expression information for tissues
    """
    tissue = models.ForeignKey('Tissue', on_delete=models.CASCADE)
    gene = models.ForeignKey('Gene', on_delete=models.CASCADE)
    expression_level = models.FloatField()  # TPM
    expression_level_norm = models.FloatField(null=True)  # log 2 normalized

    class Meta:
        unique_together = ('tissue', 'gene')


class GeneExpressionLevel(models.Model):
    """
    Expresion information for Cancer types
    """
    expression_cancer_type = models.ForeignKey('ExpressionCancerType', on_delete=models.CASCADE)
    gene = models.ForeignKey('Gene', on_delete=models.CASCADE)
    expression_level = models.FloatField()  # Batch effect log 2 normalized
    expression_level_tpm = models.FloatField(null=True)  # TPM

    class Meta:
        unique_together = ('expression_cancer_type', 'gene')


class MutationCancerType(models.Model):
    name = models.CharField(max_length=128, default='')
    abbreviation = models.CharField(max_length=128, default='')


class MutationCounts(models.Model):
    cancer_type = models.ForeignKey('MutationCancerType', on_delete=models.CASCADE)
    gene = models.ForeignKey('Gene', on_delete=models.CASCADE)
    mutation_counts = models.FloatField()
    mutation_score = models.FloatField()

    class Meta:
        unique_together = ('cancer_type', 'gene')


class CancerGeneEntity(models.Model):
    cancer_dataset_id = models.ForeignKey('CancerDataset', on_delete=models.CASCADE)
    cancer_type_id = models.ForeignKey('CancerType', on_delete=models.CASCADE)
    cancer_gene_id = models.ForeignKey('Gene', on_delete=models.CASCADE)
    pubmed_id = models.PositiveIntegerField(default=None, null=True)
    cancer_occurrences = models.PositiveIntegerField(default=None, null=True)

    class Meta:
        unique_together = ('cancer_dataset_id', 'cancer_type_id', 'cancer_gene_id')


class MinSpanningTree(models.Model):
    cancer_dataset = models.ForeignKey('CancerDataset', on_delete=models.CASCADE)
    edge_dataset = models.ForeignKey('InteractionGeneGeneDataset', on_delete=models.CASCADE)
    cancer_type = models.ForeignKey('CancerType', on_delete=models.CASCADE)
    gene = models.ForeignKey('Gene', on_delete=models.CASCADE)
    is_cancer = models.BooleanField(default=False)

    class Meta:
        unique_together = ('cancer_dataset', 'edge_dataset', 'cancer_type', 'gene')


class DrugTargetAction(models.Model):
    name = models.CharField(max_length=128, default='', unique=True)


class GeneDrugInteraction(models.Model):
    gene_id = models.ForeignKey('Gene', on_delete=models.CASCADE, related_name='relation_gene_id')
    drug_id = models.ForeignKey('Drug', on_delete=models.CASCADE, related_name='relation_gene_drug_id')
    interaction_dataset_id = models.ForeignKey('InteractionGeneDrugDataset', on_delete=models.CASCADE, unique=False)
    pubmed_id = models.CharField(max_length=256, default='', null=True)
    action = models.CharField(max_length=50, default='', null=True)

    def __str__(self) -> str:
        return f'{self.interaction_dataset_id.name}-{self.drug_id.db_id}-{self.gene_id.name}'

    class Meta:
        unique_together = ('gene_id', 'drug_id', 'interaction_dataset_id')


class RelationDatasetType(models.Model):
    cancer_dataset_id = models.ForeignKey('CancerDataset', on_delete=models.CASCADE)
    cancer_type_id = models.ForeignKey('CancerType', on_delete=models.CASCADE)

    class Meta:
        unique_together = ('cancer_dataset_id', 'cancer_type_id')


class InteractionGeneGeneDataset(models.Model):
    """
    Class for interaction datasets like "StringDB" or "Biogrid"

    All interaction entries have an ID to enable filtering by interaction dataset
    """
    name = models.CharField(max_length=128, default='', unique=True)
    link = models.CharField(max_length=128, default='', unique=True)
    version = models.CharField(max_length=64, default='')
    n_interactions = models.PositiveIntegerField(default=None, null=True)


class InteractionGeneDrugDataset(models.Model):
    """
    Class for interaction datasets like "StringDB" or "Biogrid"

    All interaction entries have an ID to enable filtering by interaction dataset
    """
    name = models.CharField(max_length=128, default='', unique=True)
    link = models.CharField(max_length=128, default='', unique=True)
    version = models.CharField(max_length=64, default='')
    n_interactions = models.PositiveIntegerField(default=None, null=True)


class GeneGeneInteraction(models.Model):
    gene_a = models.ForeignKey('Gene', on_delete=models.CASCADE, related_name='interaction_gene_a')
    gene_b = models.ForeignKey('Gene', on_delete=models.CASCADE, related_name='interaction_gene_b')
    interaction_dataset_id = models.ForeignKey('InteractionGeneGeneDataset', on_delete=models.CASCADE)
    pubmed_id = models.CharField(max_length=10000, default='', null=True)
    confidence_value = models.CharField(max_length=1500, default='', null=True)
    detection_method_psi_mi = models.CharField(max_length=5000, default='', null=True)
    detection_method_name = models.CharField(max_length=15000, default='', null=True)
    type_psi_mi = models.CharField(max_length=5000, default='', null=True)
    type_name = models.CharField(max_length=15000, default='', null=True)

    class Meta:
        unique_together = ('gene_a', 'gene_b', 'interaction_dataset_id')


class Disease(models.Model):
    mondo_id = models.PositiveIntegerField(unique=True)
    icd_10 = models.CharField(max_length=10, default='', null=True)
    name = models.CharField(max_length=200, default='', unique=True)

    class Meta:
        unique_together = ('mondo_id', 'icd_10')


class ShortestDistanceGeneToCancerGene(models.Model):
    cancer_dataset = models.ForeignKey('CancerDataset', on_delete=models.CASCADE)
    gene_interaction_dataset = models.ForeignKey('InteractionGeneGeneDataset', on_delete=models.CASCADE)
    cancer_type = models.ForeignKey('CancerType', on_delete=models.CASCADE, default=None)
    gene_a = models.ForeignKey('Gene', on_delete=models.CASCADE, related_name='short_dist_gene_a')
    gene_b = models.ForeignKey('Gene', on_delete=models.CASCADE, related_name='short_dist_gene_b')
    distance = models.IntegerField(default=None, blank=True, null=True)


class ShortestDistanceDrugToCancerGene(models.Model):
    cancer_dataset = models.ForeignKey('CancerDataset', on_delete=models.CASCADE)
    gene_drug_interaction_dataset = models.ForeignKey('InteractionGeneDrugDataset', on_delete=models.CASCADE)
    cancer_type = models.ForeignKey('CancerType', on_delete=models.CASCADE, default=None)
    drug = models.ForeignKey('Drug', on_delete=models.CASCADE)
    gene = models.ForeignKey('Gene', on_delete=models.CASCADE)
    distance = models.IntegerField(default=None, blank=True, null=True)


class DiseaseInteractions(models.Model):
    disease_a = models.ForeignKey('Disease', on_delete=models.CASCADE, related_name='disease_obj_a')
    disease_b = models.ForeignKey('Disease', on_delete=models.CASCADE, related_name='disease_obj_b')
    phi_cor = models.FloatField(default=None, blank=True, null=True)
    relative_risk = models.FloatField(default=None, blank=True, null=True)
    relative_risk_norm = models.FloatField(default=None, blank=True, null=True)

    class Meta:
        unique_together = ('disease_a', 'disease_b')


class DiseaseGeneInteractions(models.Model):
    disease = models.ForeignKey('Disease', on_delete=models.CASCADE)
    gene = models.ForeignKey('Gene', on_delete=models.CASCADE)
    database_assertedBy = models.CharField(max_length=128, default=None, blank=True, null=True)
    score_DisGeNet = models.FloatField(default=None, blank=True, null=True)

    class Meta:
        unique_together = ('disease', 'gene')


class EmpiricalNodeData(models.Model):
    content_type = models.ForeignKey(ContentType, on_delete=models.CASCADE)
    object_id = models.PositiveIntegerField()
    content_object = GenericForeignKey('content_type', 'object_id')
    gene_interaction_dataset = models.ForeignKey('InteractionGeneGeneDataset', on_delete=models.CASCADE)
    drug_interaction_dataset = models.ForeignKey('InteractionGeneDrugDataset', on_delete=models.CASCADE, null=True)
    input_size = models.PositiveIntegerField()
    count = models.PositiveIntegerField()

    class Meta:
        unique_together = ('object_id', 'gene_interaction_dataset', 'drug_interaction_dataset', 'input_size')


class Task(models.Model):
    token = models.CharField(max_length=32, unique=True)
    created_at = models.DateTimeField(auto_now_add=True)
    target = models.CharField(max_length=32, choices=[('drug', 'Drug'), ('drug-target', 'Drug Target')])

    algorithm = models.CharField(max_length=128)
    parameters = models.TextField()

    progress = models.FloatField(default=0.0)  # Progress as fraction (0.0 - 1.0)
    started_at = models.DateTimeField(null=True)
    finished_at = models.DateTimeField(null=True)
    worker_id = models.CharField(max_length=128, null=True)
    job_id = models.CharField(max_length=128, null=True)
    done = models.BooleanField(default=False)
    failed = models.BooleanField(default=False)
    status = models.CharField(max_length=255, null=True)

    result = models.TextField(null=True)
