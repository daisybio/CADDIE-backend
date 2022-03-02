from caddie.management.includes.DataCleaner import DataCleaner
from caddie.management.includes.DatabaseController import DatabaseController
from caddie import models


class DatabasePopulator:
    """
    Class that provides functions to populate the database.
    Uses functions of the class DataCleaner to receive cleaned data, which then can be passed
    to functions of the class DatabaseController to add it to the database.

    Functions of this class are called in the script "populate_db.py" to populate the database initially
    """
    def __init__(self):
        self.data_cleaner = DataCleaner()
        self.dbc = DatabaseController()

    def populate_gene_model(self):
        print('Populating Gene model ... \n')

        print('Loading Common Genes \n')
        df = self.data_cleaner.process_common_gene_data()
        self.dbc.add_common_gene_dataset(df)

        print('Done! \n')

    def populate_cancer_gene_model(self):
        print('Populating Cancer Gene models ...')

        # load the NCG6 dataset
        print('Loading NCG version 6 dataset')
        df = self.data_cleaner.process_ncg6_data()
        # add the complete NCG6 dataset
        self.dbc.add_gene_dataset(
            dataset_name='NCG',
            link='http://ncg.kcl.ac.uk/',
            df=df,
            version='6'
        )

        # load the NCG6 dataset
        print('Loading NCG version 7 dataset')
        df = self.data_cleaner.process_ncg7_data()
        # add the complete NCG6 dataset
        self.dbc.add_gene_dataset(
            dataset_name='NCG',
            link='http://ncg.kcl.ac.uk/',
            df=df,
            version='7'
        )

        # load COSMIC dataset
        print('Loading COSMIC dataset')
        df = self.data_cleaner.process_cosmic_data()
        # add the complete COSMIC dataset
        self.dbc.add_gene_dataset(
            dataset_name='COSMIC',
            link='https://cancer.sanger.ac.uk/cosmic',
            df=df,
            version='v92'
        )

        print('Loading IntOGen dataset')
        df = self.data_cleaner.process_intogen_cancer_genes()
        self.dbc.add_gene_dataset(
            dataset_name='IntOGen',
            link='https://www.intogen.org/search',
            df=df,
            version='2020-02-01'
        )

        print('Loading cancer-genes.org dataset')
        df = self.data_cleaner.process_cancer_genes_org_cancer_genes()
        self.dbc.add_gene_dataset(
            dataset_name='cancer-genes.org',
            link='http://www.cancer-genes.org/',
            df=df,
            version='2020-11-20'
        )

        print('Done!\n')

    def populate_drug_model(self):
        print('Populating Drug model ...')

        print('Loading Common DRUG dataset')
        # get the cleaned dataset
        df = self.data_cleaner.process_common_drug_data()
        # add it to the database
        self.dbc.add_complete_drug_dataset(
            dataset_name='Common',
            df=df)

        print('Done!\n')

    def populate_gene_gene_interactions(self):
        print('Populating Gene Gene Interaction model ...')

        print('Loading BioGRID Dataset')
        df = self.data_cleaner.process_common_protein_protein_data()
        self.dbc.add_gene_gene_relations(
            df=df,
            interaction_dataset_name='BioGRID',
            link='https://thebiogrid.org/',
            version='4.0'
        )

        # TODO remove this (files + functions)
        # print('Loading BIOGRID-ORGANISM-Homo_sapiens dataset')
        # get the cleaned dataset
        # df = self.data_cleaner.process_biogrid_target_target_data()
        # add it to the database
        # self.dbc.add_gene_gene_relations(
        #     df=df,
        #     interaction_dataset_name='BioGRID',
        #     link='https://thebiogrid.org/',
        #     version='4.0'
        # )

        print('Loading STRING Dataset')
        df = self.data_cleaner.process_stringdb_gene_gene_interactions()
        self.dbc.add_gene_gene_relations(
            df=df,
            interaction_dataset_name='STRING',
            link='https://string-db.org/',
            version='11.0'
        )

        print('Loading APID Dataset')
        df = self.data_cleaner.process_apid_interaction_data()
        self.dbc.add_gene_gene_relations(
            df=df,
            interaction_dataset_name='APID',
            link='http://cicblade.dep.usal.es:8080/APID/',
            version='January 2019'
        )

        print('Loading IID Dataset')
        df = self.data_cleaner.process_iid_interaction_data()
        self.dbc.add_gene_gene_relations(
            df=df,
            interaction_dataset_name='IID',
            link='http://iid.ophid.utoronto.ca/',
            version='2018-11'
        )

        print('Loading HTRI Dataset')
        df = self.data_cleaner.process_htri_interactions()
        self.dbc.add_gene_gene_relations(
            df=df,
            interaction_dataset_name='HTRIdb',
            link='http://www.lbbc.ibb.unesp.br/htri/index.jsp',
            version='2020-11-20'
        )

        print('Loading REACTOME Dataset')
        df = self.data_cleaner.process_reactome_interactions()
        self.dbc.add_gene_gene_relations(
            df=df,
            interaction_dataset_name='REACTOME',
            link='https://reactome.org/',
            version='74 Panther September 2020'
        )

        print('Done!\n')

    def populate_drug_gene_interactions(self):
        print('Populating Drug Gene Interaction model ...')

        # print('Loading Common Drug Protein dataset\n')
        # # get the cleaned dataset
        # df = self.data_cleaner.process_common_drug_protein_data()
        # # add it to the database
        # self.dbc.add_drug_gene_relations(
        #     df=df,
        #     interaction_dataset_name='CoVex',
        #     link='https://exbio.wzw.tum.de/covex/',
        #     version='1.0',
        #     common=True
        # )

        print('Loading BIOGRID-ORGANISM-Homo_sapiens dataset')
        # get the cleaned dataset
        df_interactions = self.data_cleaner.process_common_drug_protein_data()
        df_proteins_genes = self.data_cleaner.process_common_gene_data()
        df = df_interactions.merge(df_proteins_genes[['protein_ac', 'entrez_id']], how='inner', on='protein_ac')
        # add it to the database
        self.dbc.add_drug_gene_relations(
            df=df,
            interaction_dataset_name='BioGRID',
            link='https://thebiogrid.org/',
            version='4.0',
            common=True
        )

        print('Loading drugbank dataset')
        df = self.data_cleaner.process_drugbank_drug_gene_interactions()
        self.dbc.add_drug_gene_relations_v2(
            df=df,
            interaction_dataset_name='DrugBank',
            link='https://go.drugbank.com/',
            version='5.1.7'
        )

        print('Loading drugbank actions dataset')
        df = self.data_cleaner.process_drugbank_drug_gene_actions()
        self.dbc.add_drug_gene_actions(df, 'DrugBank')

        print('Loading ChEMBL dataset')
        df = self.data_cleaner.process_chembl_drug_gene_interactions()
        self.dbc.add_drug_gene_relations(
            df=df,
            interaction_dataset_name='ChEMBL',
            link='https://www.ebi.ac.uk/chembl/',
            version='27',
            common=True
        )

        print('Loading DGIdb dataset')
        df = self.data_cleaner.process_dgidb_drug_gene_interactions()
        self.dbc.add_drug_gene_relations_v2(
            df=df,
            interaction_dataset_name='DGIdb',
            link='https://www.dgidb.org/',
            version='4.2.0'
        )

        print('Done!\n')

    def populate_tissue_model(self):
        print('Populating Tissue and ExpressionLevel model ...')
        df = self.data_cleaner.process_tissue_data()

        self.dbc.add_tissues_and_expression(df)

        print('Normalizing Tissue and ExpressionLevel model ...')
        self.dbc.normalize_tissue_expression()

        print('Done!\n')

    def populate_exp_model(self):
        print('Populating Expression Cancer Type and Gene Expression Level model ...')
        df_expr_norm, df_expr_tpm = self.data_cleaner.process_expression_data()

        self.dbc.add_expression(df_expr_norm, df_expr_tpm)

        print('Done!\n')

    def add_cancer_type_occurrences(self):

        print('Populating cancer_occurances for genes ...\n')

        self.dbc.add_cancer_type_occurrences()

        print('Done!\n')

    def add_comorbidities(self):

        print('Populating comorbidities ...\n')

        # TODO remove this? Not used right now, maybe later
        # df = self.data_cleaner.process_comorbidity_diseases()
        # self.dbc.add_comorbidity_disease_data(df)
        #
        # df = self.data_cleaner.process_comorbidity_diseases_interactions()
        # self.dbc.add_comorbidity_disease_interactions(df)

        df = self.data_cleaner.process_comorbidity_gene_interactions()
        self.dbc.add_comorbidity_gene_interactions(df)

        print('Done!\n')

    def add_shortest_distances(self):

        print('Populating shortest distances ...\n')

        print('Shortest distances for Cancer Genes\n')

        cancer_datasets = [e.name for e in models.CancerDataset.objects.all()]
        gene_interaction_datasets = [e.name for e in models.InteractionGeneGeneDataset.objects.all()]
        drug_gene_interaction_datasets = [e.name for e in models.InteractionGeneDrugDataset.objects.all()]

        for cancer_dataset in cancer_datasets:
            for gene_interaction_dataset in gene_interaction_datasets:
                for drug_gene_interaction_dataset in drug_gene_interaction_datasets:
                    df = self.data_cleaner.process_shortest_distances(
                        cancer_dataset,
                        gene_interaction_dataset,
                        drug_gene_interaction_dataset
                    )

                    self.dbc.add_shortest_distances(
                        df,
                        cancer_dataset,
                        gene_interaction_dataset,
                        drug_gene_interaction_dataset
                    )

        print('Done!\n')

    def add_mutation_count(self):

        print('Populating mutation counts ...\n')

        df = self.data_cleaner.process_mutation_counts()
        self.dbc.add_mutations_counts(df)

        print('Done!\n')

    def add_cancernet(self):

        print('Populating cancernet counts ...\n')
        regions = ['NIH', 'UK']
        for region in regions:
            df = self.data_cleaner.process_cancernet_targeted(region)
            self.dbc.add_cancernet_table(df, region, True)
        for region in regions:
            df = self.data_cleaner.process_cancernet_untargeted(region)
            self.dbc.add_cancernet_table(df, region, False)

        print('Done!\n')
