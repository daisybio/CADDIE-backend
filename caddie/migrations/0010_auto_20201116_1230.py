# Generated by Django 3.0.5 on 2020-11-16 12:30

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('caddie', '0009_gene_n_mutations'),
    ]

    operations = [
        migrations.AlterField(
            model_name='disease',
            name='mondo_id',
            field=models.PositiveIntegerField(unique=True),
        ),
    ]
