# Generated by Django 4.2.4 on 2024-04-26 10:36

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0011_alter_sample_well'),
    ]

    operations = [
        migrations.AlterField(
            model_name='sample',
            name='sample_number',
            field=models.TextField(blank=True, null=True),
        ),
    ]
