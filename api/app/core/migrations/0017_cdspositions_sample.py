# Generated by Django 4.2.4 on 2024-07-16 07:35

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0016_remove_samplecount_consensus_fraction_n_cds_and_more'),
    ]

    operations = [
        migrations.AddField(
            model_name='cdspositions',
            name='sample',
            field=models.ForeignKey(default=None, on_delete=django.db.models.deletion.CASCADE, to='core.sample'),
            preserve_default=False,
        ),
    ]
