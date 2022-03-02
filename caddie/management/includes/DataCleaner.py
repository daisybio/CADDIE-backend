import glob
import os
import pandas as pd
import numpy as np
import urllib.parse
import urllib.request
from io import StringIO
import html


class DataCleaner:
    """
    This class provides functions to clean data.
    Functions are mostly customized to specific datasets.
    Functions:
        1. read data files
        2. clean the data
        3. return cleaned data in a format that can be picked up by DatabasePopulator to add it to the db

    Made for the use in populate_db to initialize the database in the beginning.
    """

    def read_file(self, filename, names=None, dtype=None, index_col=None, sep=None):
        """
        dynamically search the data folder to load specific files
        """
        data_folders = glob.glob('data/*')
        # iterating over all folders
        for data_folder in data_folders:
            # iterating over all files
            for file in os.listdir(data_folder):
                name = file.split('.')[0]
                # file is file we are searching for
                if name == filename:
                    file_path = os.path.join(data_folder, file)
                    file_path_abs = os.path.abspath(file_path)

                    if sep is None:
                        # comma sep is default
                        sep = ','

                        # all tab separated files
                        if file[-3:] == 'tsv' or file[-7:] == 'tab.txt' or file[-3:] == 'tab':
                            sep = '\t'
                    return pd.read_csv(
                        file_path_abs, sep=sep, names=names, dtype=dtype,
                        index_col=index_col, compression='gzip', header=0)

    def process_ncg6_data(self):
        """
        Reads and cleans the NCG6 dataset.
        The NCG6 dataset contains cancer-driver genes with attributes:
        - entrez
        - symbol
        - pubmed_id
        - type
        - primary_site
        - cancer_type
        - method

        :return: pandas DataFrame with the relevant data ('name', 'cancer_type', 'entrez_id', 'pubmed_id')
        """

        df = self.read_file('NCG6_cancergenes')[['symbol', 'cancer_type', 'entrez', 'pubmed_id']]

        df = df.drop_duplicates()

        # data is very unclean, but wtf is this cancer type... No cancer???
        df = df[df['cancer_type'] != '-']

        # convert html entities like &#39; to string
        df['cancer_type'] = df['cancer_type'].map(html.unescape)

        df = df.rename(columns={
            'symbol': 'name',
            'entrez': 'entrez_id',
        })

        df = df.drop_duplicates()

        return df

    def process_ncg7_data(self):
        """
        Reads and cleans the NCG7 dataset.
        The NCG6 dataset contains cancer-driver genes with attributes:
        - entrez
        - symbol
        - pubmed_id
        - type
        - primary_site
        - cancer_type
        - method

        :return: pandas DataFrame with the relevant data ('name', 'cancer_type', 'entrez_id', 'pubmed_id')
        """

        df = self.read_file('NCG7_cancergenes')[['symbol', 'cancer_type', 'entrez', 'pubmed_id']]

        df = df.drop_duplicates()

        # data is very unclean, but wtf is this cancer type... No cancer???
        df = df[df['cancer_type'] != '-']

        # convert html entities like &#39; to string
        df['cancer_type'] = df['cancer_type'].map(html.unescape)

        df = df.rename(columns={
            'symbol': 'name',
            'entrez': 'entrez_id',
        })

        df = df.drop_duplicates()

        return df

    def process_cosmic_data(self):
        """
        The COSMIC Dataset is really dirty, entrez ids are not correct, more than 400 different cancer types

        Reads and cleans the COSMIC dataset.
        The COSMIC dataset contains cancer-driver genes with arributes:
        'Gene Symbol', 'Name', 'Entrez GeneId', 'Genome Location', 'Tier',
        'Hallmark', 'Chr Band', 'Somatic', 'Germline', 'Tumour Types(Somatic)',
        'Tumour Types(Germline)', 'Cancer Syndrome', 'Tissue Type',
        'Molecular Genetics', 'Role in Cancer', 'Mutation Types',
        'Translocation Partner', 'Other Germline Mut', 'Other Syndrome',
        'Synonyms'
        :return: pandas DataFrame with the relevant data ['name', 'entrez_id', 'cancer_type']]

        """

        df = self.read_file('cancer_gene_census')
        df['Tumour Types'] = df['Tumour Types(Somatic)'].astype(str) + ', ' + df['Tumour Types(Germline)'].astype(str)

        # convert 'Tumour Types(Somatic)' into list
        df['Tumour Types'] = df['Tumour Types'].astype(str)
        df['Tumour Types'] = df['Tumour Types'].apply(lambda x: x.split(', '))

        df = df.explode('Tumour Types')
        df['Tumour Types'] = df['Tumour Types'].astype(str).replace('nan', np.nan)

        df = df.dropna(subset=['Tumour Types'])
        df['Tumour Types'] = df['Tumour Types'].apply(lambda x: x.replace(' ', '_'))

        # manually fill in entrez ids for the 4 genes with missing ids
        fill_in_gene_ids = [('MALAT1', 378938), ('DUX4L1', 22947), ('HMGN2P46', 283651), ('MDS2', 259283)]
        for gene, entrez_id in fill_in_gene_ids:
            index = df[df['Gene Symbol'] == gene].index[0]
            df.at[index, 'Entrez GeneId'] = entrez_id

        # data is very unclean, but wtf is this cancer type... No cancer???
        df = df[df['Tumour Types'] != '-']

        # no pubmed_id
        df['pubmed_id'] = None

        df = df[['Gene Symbol', 'Entrez GeneId', 'Tumour Types', 'pubmed_id']]

        # cosmic dataset has dirty entrez ids, we need to fill in missing entrez ids and replace wrong ones
        fill_in_gene_ids = [
            ('MALAT1', 378938),
            ('DUX4L1', 22947),
            ('HMGN2P46', 283651),
            ('MDS2', 259283),
            ('EIF1AX', 1964)]

        # fill in the entrez ids
        for gene, entrez_id in fill_in_gene_ids:
            index = df[df['Gene Symbol'] == gene].index[0]
            df.at[index, 'Entrez GeneId'] = entrez_id

        df = df.rename(columns={
            'Gene Symbol': 'name',
            'Entrez GeneId': 'entrez_id',
            'Tumour Types': 'cancer_type'
        })

        df = df.drop_duplicates()

        return df

    def process_common_drug_data(self):
        """
        :return: pandas DataFrame with columns, drug_id, drug_name, drug_status, in_trial, in_literature, links
        """

        df = self.read_file('drug-intrial-link-file')

        # dont use the CoVex status, use new status from repo trial which is not just Rona-specific
        df = df.drop(['drug_status'], axis=1)

        # merge 'is_nutraceutical'
        df_drugclasses_repodb = self.read_file('drugbank_drug_groups')
        df_drugclasses_repodb['drug_id'] = df_drugclasses_repodb['primaryDomainId'].map(lambda x: x.split('.')[1])

        df = df.merge(df_drugclasses_repodb[['drug_id', 'nutraceutical', 'status']], on='drug_id', how='left')

        df['nutraceutical'] = df['nutraceutical'].fillna(False)

        # read atc drug list and add boolean if drug is in list
        df_atc = self.read_file('atc_antineoplastic_and_immunomodulating_agents')
        df['is_atc_antineoplastic_and_immunomodulating_agent'] = df['drug_id'].isin(df_atc['id'])

        df['is_atc_antineoplastic_and_immunomodulating_agent'] = \
            df['is_atc_antineoplastic_and_immunomodulating_agent'].fillna(False)

        #  adding pharmaGKD_id 	pubchem_substance_id
        df_identifier = self.read_file('drugbank_labels')

        df = df.merge(df_identifier, on='drug_id', how='left')

        df = df[df['status'].notna()]

        return df

    def process_stringdb_gene_gene_interactions(self):
        """

        :return: pandas DataFrame with columns: 'protein1', 'protein2', 'neighborhood', 'neighborhood_transferred',
       'fusion', 'cooccurence', 'homology', 'coexpression',
       'coexpression_transferred', 'experiments', 'experiments_transferred',
       'database', 'database_transferred', 'textmining',
       'textmining_transferred', 'combined_score', 'entrez_a', 'entrez_b'
        """

        df = self.read_file('string_interactions')

        return df

    def process_biogrid_drug_data(self):
        """
        Reads and cleans the BIOGRID CHEMICALS dataset.

        :return: pandas DataFrame with the relevant data
        """

        df = self.read_file('BIOGRID-CHEMICALS-3').drop_duplicates(subset=['Chemical Name'])

        df = df[
            [
                'Chemical Name',
                'Action',
                'Chemical Synonyms',
                'Molecular Formula',
                'Chemical Type'
            ]
        ]

        df = df.rename(columns={'Chemical Name': 'name',
                                'Action': 'action',
                                'Chemical Synonyms': 'synonym',
                                'Molecular Formula': 'molecular_formula',
                                'Chemical Type': 'type'})

        # standardized None value
        df = df.replace('unknown', None)

        return df

    def process_biogrid_target_target_data(self):
        """
        :return: pandas Dataframe with columns: 'entrez_a','entrez_b','pubmed_id',
        confidence_value','detection_method_psi_mi','detection_method_name','type_psi_mi', 'type_name'
        """

        df = self.read_file('BIOGRID-ORGANISM-Homo_sapiens-3')[
            ['#ID Interactor A', 'ID Interactor B', 'Interaction Detection Method', 'Publication Identifiers',
             'Interaction Types', 'Confidence Values']]

        def parse_interactor(x):
            # format entrez protein/locuslink:6416
            # wanted: 6416
            return x.split(':')[-1]

        def parse_publication_identifiers(x):
            # format: pubmed:9006895
            # wanted: 9006895
            return x.split(':')[-1]

        def parse_interaction_detection_method_get_id(x):
            # format: psi-mi:"MI:0018"(two hybrid)
            # wanted: MI:0018
            return x.split('"')[1]

        def parse_interaction_detection_method_get_name(x):
            # format: psi-mi:"MI:0018"(two hybrid)
            # wanted: two hybrid
            return x.split('(')[1][:-1]

        def parse_interaction_types_get_id(x):
            # format: psi-mi:"MI:0407"(direct interaction)
            # wanted: MI:0407
            return x.split('"')[1]

        def parse_interaction_types_get_name(x):
            # format: psi-mi:"MI:0407"(direct interaction)
            # wanted: direct interaction
            return x.split('(')[1][:-1]

        def parse_confidence_value(x):
            # format: score:7.732982515 or '-'
            # wanted: 7.732982515 or '-'
            if x == '-':
                return '-'
            else:
                return x.split(':')[1]

        df['#ID Interactor A'] = df['#ID Interactor A'].map(parse_interactor)
        df['ID Interactor B'] = df['ID Interactor B'].map(parse_interactor)

        df['Interaction Detection Method ID'] = df['Interaction Detection Method'].map(
            parse_interaction_detection_method_get_id)
        df['Interaction Detection Method Name'] = df['Interaction Detection Method'].map(
            parse_interaction_detection_method_get_name)
        df = df.drop('Interaction Detection Method', axis=1)

        df['Publication Identifiers'] = df['Publication Identifiers'].map(parse_publication_identifiers)

        df['Interaction Types ID'] = df['Interaction Types'].map(parse_interaction_types_get_id)
        df['Interaction Types Name'] = df['Interaction Types'].map(parse_interaction_types_get_name)
        df = df.drop('Interaction Types', axis=1)

        df['Confidence Values'] = df['Confidence Values'].map(parse_confidence_value)

        # remove dirty data (entrez id is sometimes protein name
        to_remove = ['P0DTC1', 'P0DTD2', 'Q7TLC7']
        df = df[~df['#ID Interactor A'].isin(to_remove)]

        df = df.rename(
            columns={
                '#ID Interactor A': 'entrez_a',
                'ID Interactor B': 'entrez_b',
                'Publication Identifiers': 'pubmed_id',
                'Confidence Values': 'confidence_value',
                'Interaction Detection Method ID': 'detection_method_psi_mi',
                'Interaction Detection Method Name': 'detection_method_name',
                'Interaction Types ID': 'type_psi_mi',
                'Interaction Types Name': 'type_name'
            }
        )

        return df

    def process_biogrid_drug_target_data(self):
        """
        reads BIOGRID-CHEMICALS-3.5.187.chemtab.txt file and parses it into formatted output
        so that i can
        :return: pandas Dataframe with column names 'Chemical Name', 'Entrez Gene ID', 'Pubmed ID'
        """
        df = self.read_file('BIOGRID-CHEMICALS-3')[['Chemical Name', 'Entrez Gene ID', 'Pubmed ID']]

        df = df.drop_duplicates()

        df = df.rename(columns={'Chemical Name': 'drug_name',
                                'Entrez Gene ID': 'entrez_id',
                                'Pubmed ID': 'pubmed_id'
                                })

        return df

    def process_biogrid_gene_data(self):

        df = self.read_file('BIOGRID-ORGANISM-Homo_sapiens-3')[
            ['#ID Interactor A', 'ID Interactor B', 'Alt IDs Interactor A', 'Alt IDs Interactor B',
             'Aliases Interactor A', 'Aliases Interactor B']]

        def parse_interactor(x):
            # format entrez protein/locuslink:6416
            # wanted: 6416
            return x.split(':')[-1]

        def parse_gene_name(x):
            # 'format: biogrid:112315|entrez protein/locuslink:MAP2K4
            # |uniprot/swiss-prot:P45985|refseq:NP_003001|refseq:NP_001268364'

            # wanted: MAP2K4

            # x1 = entrez protein/locuslink:MAP2K4
            x1 = x.split('|')[1]

            return x1.split(':')[1]

        def parse_aliase_interactors(x):
            # format: entrez protein/locuslink:JNKK(protein name synonym)|
            # entrez protein/locuslink:JNKK1(protein name synonym)|
            # entrez protein/locuslink:MAPKK4(protein name synonym)|
            # entrez protein/locuslink:MEK4(protein name synonym)|
            # entrez protein/locuslink:MKK4(protein name synonym)|
            # entrez protein/locuslink:PRKMK4(protein name synonym)|
            # entrez protein/locuslink:SAPKK-1(protein name synonym)|
            # entrez protein/locuslink:SAPKK1(protein name synonym)|
            # entrez protein/locuslink:SEK1(protein name synonym)|
            # entrez protein/locuslink:SERK1(protein name synonym)|
            # entrez protein/locuslink:SKK1(protein name synonym)

            # wanted: JNKK,JNKK1,MAPKK4,...

            x_cut1 = x.replace('entrez protein/locuslink:', '')
            x_cut2 = x_cut1.replace('(protein name synonym)', '')
            x_cut3 = x_cut2.replace(' ', '')
            x_cut4 = x_cut3.replace('-', '')

            names = x_cut4.split('|')

            return ','.join(names)

        df['#ID Interactor A'] = df['#ID Interactor A'].map(parse_interactor)
        df['ID Interactor B'] = df['ID Interactor B'].map(parse_interactor)

        df['Alt IDs Interactor A'] = df['Alt IDs Interactor A'].map(parse_gene_name)
        df['Alt IDs Interactor B'] = df['Alt IDs Interactor B'].map(parse_gene_name)

        df['Aliases Interactor A'] = df['Aliases Interactor A'].map(parse_aliase_interactors)
        df['Aliases Interactor B'] = df['Aliases Interactor B'].map(parse_aliase_interactors)

        df = df.rename(
            columns={
                '#ID Interactor A': 'entrez_a',
                'ID Interactor B': 'entrez_b',
                'Alt IDs Interactor A': 'gene_name_a',
                'Alt IDs Interactor B': 'gene_name_b',
                'Aliases Interactor A': 'alias_gene_a',
                'Aliases Interactor B': 'alias_gene_b'
             }
        )

        df_a = df[['entrez_a', 'gene_name_a', 'alias_gene_a']]
        df_b = df[['entrez_b', 'gene_name_b', 'alias_gene_b']]

        df_a = df_a.rename(
            columns={
                'entrez_a': 'entrez_id',
                'gene_name_a': 'name',
                'alias_gene_a': 'alias'
            }
        )

        df_b = df_b.rename(
            columns={
                'entrez_b': 'entrez_id',
                'gene_name_b': 'name',
                'alias_gene_b': 'alias'
            }
        )

        df = pd.concat([df_a, df_b])

        # remove all dirty entries where entrez_id is not a real entrez_id
        # entrez_ids are numbers
        df = df[df['entrez_id'].str.isnumeric()]

        df = df.drop_duplicates('name')

        return df

    def entrez_to_uniprot(self, entrez_ids: list):
        # make sure entrez_ids are strings

        entrez_ids = list(map(str, entrez_ids))

        columns = [
            'UniProtKB-AC',
            'UniProtKB-ID',
            'GeneID (EntrezGene)',
            'RefSeq',
            'GI',
            'PDB',
            'GO',
            'UniRef100',
            'UniRef90',
            'UniRef50',
            'UniParc',
            'PIR',
            'NCBI-taxon',
            'MIM',
            'UniGene',
            'PubMed',
            'EMBL',
            'EMBL-CDS',
            'Ensembl',
            'Ensembl_TRS',
            'Ensembl_PRO',
            'Additional PubMed'
        ]

        df = self.read_file('HUMAN_9606_idmapping_selected', names=columns, dtype='str')

        # get sub_df containing only rows with given entrez_ids
        df_sub = df[df['GeneID (EntrezGene)'].isin(entrez_ids)]

        # create a lookup series with entrez_id: uniprot_kb_ac relations
        lookup_series = df_sub[['UniProtKB-AC', 'GeneID (EntrezGene)']].set_index('GeneID (EntrezGene)')

        return lookup_series

    def uniprot_to_entrez(self, uniprot_ids: list):
        # make sure uniprot_ids are strings

        uniprot_ids = list(map(str, uniprot_ids))

        columns = [
            'UniProtKB-AC',
            'UniProtKB-ID',
            'GeneID (EntrezGene)',
            'RefSeq',
            'GI',
            'PDB',
            'GO',
            'UniRef100',
            'UniRef90',
            'UniRef50',
            'UniParc',
            'PIR',
            'NCBI-taxon',
            'MIM',
            'UniGene',
            'PubMed',
            'EMBL',
            'EMBL-CDS',
            'Ensembl',
            'Ensembl_TRS',
            'Ensembl_PRO',
            'Additional PubMed'
        ]

        df = self.read_file('HUMAN_9606_idmapping_selected', names=columns, dtype='str')

        # get sub_df containing only rows with given uniprot_ids
        df_sub = df[df['UniProtKB-AC'].isin(uniprot_ids)]

        # create a lookup series with uniprot_id: entrez_id relations
        lookup_series = df_sub[['UniProtKB-AC', 'GeneID (EntrezGene)']].set_index('UniProtKB-AC')

        return lookup_series

    def uniprot_to_protein_names(self, uniprot_kb_list: list):
        # fetch all reviewed, human related entries from uniprot in tab format
        url = 'https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:9606&format=tab'
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req) as f:
            response = f.read()
        data = response.decode('utf-8')

        data_string = StringIO(data)

        df = pd.read_csv(data_string, sep='\t', dtype='str')
        # columns are: Entry, Entry name, Status, Protein names, Gene names, Organism, Length

        # filter entries to only rows with uniprot_kb:acs from uniprot_kb_list
        df_sub = df[df['Entry'].isin(uniprot_kb_list)]

        # create lookup_series with uniprot_kb_ac: Protein_names relations
        lookup_series = df_sub[['Entry', 'Protein names']].set_index('Entry')

        return lookup_series

    def entrez_to_protein(self, entrez_ids: list):
        """
        This function uses a pipline consisting of entrez_to_uniprot and uniprot_to_protein_names
        to retrieve proteins for entrez_ids
        :param entrez_ids:
        :return:
        """

        # 1 read HUMAN_9606_idmapping_selected.tab to map entrez id to UniProtKB-ID
        entrez_to_uniprot_series = self.entrez_to_uniprot(entrez_ids)

        uniprot_ids = entrez_to_uniprot_series['UniProtKB-AC'].to_list()

        # 2 lookup uniprot_id to find protein names in data downloaded from uniprot
        uniprot_to_protein_series = self.uniprot_to_protein_names(uniprot_ids)

        # 3 check fo reach entrez id if uniprotkb exists and if so, if protein names exist, else return None

        # relations stores [entrez, uniprot_acs, protein_name] lists
        relations = []
        for entrez_id in entrez_ids:

            # make sure entrez_id is string
            entrez_id = str(entrez_id)

            # try to fetch uniprot id
            try:
                uniprot_acs = entrez_to_uniprot_series.loc[entrez_id]['UniProtKB-AC'].to_list()

                for uniprot_ac in uniprot_acs:

                    # try to fetch protein names
                    try:
                        protein_name = uniprot_to_protein_series.loc[uniprot_ac]['Protein names']
                        relations.append([entrez_id, uniprot_ac, protein_name])
                    except KeyError:
                        relations.append([entrez_id, uniprot_ac, None])
                        pass

            except KeyError:
                relations.append([entrez_id, None, None])

        return pd.DataFrame(relations, columns=['entrez_id', 'uniprot_ac', 'protein'])

    def process_common_gene_data(self):
        """
        The common dataset folder contains the data used for CoVex
        :return: pandas dataframe with columns gene_name, protein_ac, protein_name, entrez_id
        """

        df = self.read_file('protein_list')
        return df

    def process_common_protein_protein_data(self):

        df = self.read_file('protein-protein-interaction')

        return df

    def process_common_drug_protein_data(self):

        df = self.read_file('drug-protein-interaction')

        return df

    def process_comorbidity_diseases(self):

        df = self.read_file('Comorbidity_ALL_with_names_edited', index_col=0)

        # we dont have name for ~40k out of ~450k icd_10 codes
        df = df.dropna()

        disease_list_a = df[['disease1_ICD10', 'disease1_Mondo', 'disease1_name']]\
            .drop_duplicates('disease1_ICD10').values
        disease_list_b = df[['disease2_ICD10', 'disease2_Mondo', 'disease2_name']]\
            .drop_duplicates('disease2_ICD10').values

        disease_list_all = [*disease_list_a, *disease_list_b]

        df_res = pd.DataFrame(disease_list_all, columns=['ICD_10', 'Mondo', 'name'])

        return df_res

    def process_comorbidity_diseases_interactions(self):
        # 'disease1_ICD10', 'disease2_ICD10', 'phi_cor', 'relative_risk',
        #        'disease1_name', 'disease2_name', 'relative_risk_norm'
        df = self.read_file('Comorbidity_ALL_with_names_edited', index_col=0)

        return df

    def process_comorbidity_gene_interactions(self):
        df = self.read_file('comorbidity_disease_gene_interactions', index_col=0)

        # filter out disease 'cancer'.... no value
        df = df[df['Name'] != 'cancer']
        return df

    def process_shortest_distances(
            self,
            cancer_dataset,
            gene_interaction_dataset,
            drug_gene_interaction_dataset
    ):
        """
        :return: pandas DataFrame with columns: source, target, dist, cancer_type
        """
        file = 'internal_{cancer_dataset}_{gene_interaction_dataset}_{drug_gene_interaction_dataset}' \
               '_shortest_distances_to_cancer_gene'.format(
                cancer_dataset=cancer_dataset,
                gene_interaction_dataset=gene_interaction_dataset,
                drug_gene_interaction_dataset=drug_gene_interaction_dataset
                )

        df = self.read_file(file)

        return df

    def process_tissue_data(self):
        """
        Returns df with expression values in fields for genes in column "UniprotId"
        :return:
        """
        df = self.read_file('expression_per_tissue_normalized', index_col=0)
        df = df.rename(columns={'Description': 'UniprotId'})
        return df

    def process_expression_data(self):
        """
        Returns df with expression values in fields for genes in column "UniprotId"
        :return:
        """
        df_expr_norm = self.read_file('cancer_expression_norm', index_col=0)
        df_expr_tpm = self.read_file('cancer_expression_tpm', index_col=0)
        return df_expr_norm, df_expr_tpm

    def process_mutation_counts(self):
        """
        loads file "gene_mutation_counts.csv"
        columns are "entrez_id" and "mutation_count"
        :return:
        """
        df = self.read_file('gene_mutation_counts')
        # not all genes have cancer types assigned
        df = df[~df['name'].isna()]
        return df

    def process_apid_interaction_data(self):
        """InteractionID	UniprotID_A	UniprotName_A	GeneName_A	UniprotID_B	UniprotName_B	GeneName_B	ExpEvidences
        Methods	Publications	3DStructures	CurationEvents"""

        df = self.read_file('9606_Q2', sep='\t', index_col=0)
        df = df.rename(columns={
            'UniprotID_A': 'from_protein_ac',
            'UniprotID_B': 'to_protein_ac'
        })
        return df[['from_protein_ac', 'to_protein_ac']]

    def process_iid_interaction_data(self):
        """
        'uniprot1',
        'uniprot2',
        'symbol1',
        'symbol2'
        ...

        NOTE we use only interactions that occur in some cancer type
        """
        df = self.read_file('human_annotated_PPIs_formatted')
        df = df.rename(columns={
            'uniprot1': 'from_protein_ac',
            'uniprot2': 'to_protein_ac'
        })
        return df

    def process_htri_interactions(self):
        df1 = self.read_file('gene-regulation-interactions')
        df1 = df1.rename(columns={'SYMBOL_TF': 'gene_name_a', 'SYMBOL_TG': 'gene_name_b'})
        df1 = df1[['gene_name_a', 'gene_name_b']]

        df2 = self.read_file('protein-protein-interactions', sep='\t',
                             names=['gene_name_a', 'gene_name_b', 'interaction'])
        df2 = df2[['gene_name_a', 'gene_name_b']]

        df = df1.append(df2).drop_duplicates()

        return df

    def process_reactome_interactions(self):

        df = self.read_file('homo-sapiens-protein-interactions', sep='\t')

        df['# Interactor 1 uniprot id'] = df['# Interactor 1 uniprot id'].map(lambda x: x.split(':')[1])
        df['Interactor 2 uniprot id'] = df['Interactor 2 uniprot id'].map(lambda x: x.split(':')[1])

        df = df.rename(
            columns={'# Interactor 1 uniprot id': 'from_protein_ac', 'Interactor 2 uniprot id': 'to_protein_ac'})

        return df

    def process_intogen_cancer_genes(self):

        df = self.read_file('intogen_cancer_genes')
        df = df.dropna()
        return df

    def process_cancer_genes_org_cancer_genes(self):
        df = self.read_file('cancer-genes-org')

        return df

    def process_drugbank_drug_gene_interactions(self):
        df = self.read_file('drugbank_drug_gene_interactions')
        df = df.dropna()
        return df

    def process_chembl_drug_gene_interactions(self):
        df = self.read_file('chembl_drug_gene_interactions')
        return df

    def process_dgidb_drug_gene_interactions(self):
        df = self.read_file('DGIdb_drug_gene_interactions')
        return df

    def process_drugbank_drug_gene_actions(self):
        df = self.read_file('drug_protein_interaction_actions')
        return df

    def process_cancernet_targeted(self, region):
        df = self.read_file(f'Drug_Target_cancertype_db-V2 - {region}')
        if 'Combination formation' not in df:
            df['Combination formation'] = ''
        return df

    def process_cancernet_untargeted(self, region):
        df = self.read_file(f'Data not included in the db-V2 - {region}')
        if 'Cancer type (from website)' not in df:
            df['Cancer type (from website)'] = df['Cancer type']
            df['Cancer type (short version)'] = ''
        return df
