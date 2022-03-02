# Generated by Django 3.0.5 on 2020-10-17 23:40

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('caddie', '0003_auto_20201017_1218'),
    ]

    operations = [
        migrations.CreateModel(
            name='ShortestDistanceDrugToCancerGene',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('distance', models.IntegerField(blank=True, default=None, null=True)),
                ('cancer_dataset', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='caddie.CancerDataset')),
                ('cancer_type', models.ForeignKey(default=None, on_delete=django.db.models.deletion.CASCADE, to='caddie.CancerType')),
                ('drug', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='caddie.Drug')),
                ('gene', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='caddie.Gene')),
                ('gene_drug_interaction_dataset', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='caddie.InteractionGeneDrugDataset')),
            ],
        ),
        migrations.CreateModel(
            name='ShortestDistanceGeneToCancerGene',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('distance', models.IntegerField(blank=True, default=None, null=True)),
                ('cancer_dataset', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='caddie.CancerDataset')),
                ('cancer_type', models.ForeignKey(default=None, on_delete=django.db.models.deletion.CASCADE, to='caddie.CancerType')),
                ('gene_a', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='short_dist_gene_a', to='caddie.Gene')),
                ('gene_b', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='short_dist_gene_b', to='caddie.Gene')),
                ('gene_interaction_dataset', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='caddie.InteractionGeneGeneDataset')),
            ],
        ),
        migrations.DeleteModel(
            name='ShortestDistanceCancerGene',
        ),
    ]
