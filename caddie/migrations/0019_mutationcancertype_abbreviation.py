# Generated by Django 3.0.5 on 2021-03-18 16:17

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('caddie', '0018_auto_20210317_2135'),
    ]

    operations = [
        migrations.AddField(
            model_name='mutationcancertype',
            name='abbreviation',
            field=models.CharField(default='', max_length=128),
        ),
    ]
