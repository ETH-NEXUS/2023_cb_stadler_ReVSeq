# Generated by Django 4.2.4 on 2024-04-26 08:36

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0009_remove_samplecount_dp1_remove_samplecount_dp2_and_more'),
    ]

    operations = [
        migrations.AddField(
            model_name='sample',
            name='control',
            field=models.BooleanField(blank=True, default=False, null=True),
        ),
        migrations.AddField(
            model_name='sample',
            name='control_type',
            field=models.CharField(blank=True, max_length=10, null=True),
        ),
    ]
