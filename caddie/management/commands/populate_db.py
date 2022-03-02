from django.core.management.base import BaseCommand
from caddie.management.includes.DatabasePopulator import DatabasePopulator
from caddie.management.includes.DatabaseController import DatabaseController


class Command(BaseCommand):
    help = 'Functions to populate the db'

    def add_arguments(self, parser):

        # dataset directory
        parser.add_argument('-d', '--drug', type=str, help='Viral protein file')
        parser.add_argument('-cg', '--add_data', type=str, help='Add data to db')
        parser.add_argument('-dm', '--clear_model', type=str, help='Delete all entries of a model')

    def handle(self, *args, **kwargs):

        db_populator = DatabasePopulator()
        dbc = DatabaseController()

        if kwargs['clear_model'] is not None:
            model_list = kwargs['clear_model'].split(',')
            dbc.delete_models(model_list)
            return

        if kwargs['add_data'] is not None:
            to_add = kwargs['add_data'].split(',')

            for data in to_add:
                if data == 'cancer_genes':
                    db_populator.populate_cancer_gene_model()
                elif data == 'drugs':
                    db_populator.populate_drug_model()
                elif data == 'gene_gene_interactions':
                    db_populator.populate_gene_gene_interactions()
                elif data == 'drug_gene_interactions':
                    db_populator.populate_drug_gene_interactions()
                elif data == 'genes':
                    db_populator.populate_gene_model()
                elif data == 'tissues':
                    db_populator.populate_tissue_model()
                elif data == 'cancer_expression':
                    db_populator.populate_exp_model()
                elif data == 'cancer_occurrences':
                    db_populator.add_cancer_type_occurrences()
                elif data == 'comorbidities':
                    db_populator.add_comorbidities()
                elif data == 'shortest_distances':
                    db_populator.add_shortest_distances()
                elif data == 'mutation_count':
                    db_populator.add_mutation_count()
                elif data == 'cancernet':
                    db_populator.add_cancernet()
