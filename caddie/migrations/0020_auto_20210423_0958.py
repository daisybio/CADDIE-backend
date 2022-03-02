# Generated by Django 3.0.5 on 2021-04-23 09:58

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('caddie', '0019_mutationcancertype_abbreviation'),
    ]

    operations = [
        migrations.CreateModel(
            name='ExpressionCancerType',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(default='', max_length=128)),
                ('abbreviation', models.CharField(default='', max_length=128)),
            ],
        ),

        migrations.CreateModel(
            name='GeneExpressionLevel',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('expression_level', models.FloatField()),
                ('expression_cancer_type', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='caddie.ExpressionCancerType')),
                ('gene', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='caddie.Gene')),
            ],
            options={
                'unique_together': {('expression_cancer_type', 'gene')},
            },
        ),
    ]
