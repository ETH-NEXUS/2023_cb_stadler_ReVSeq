# Generated by Django 4.2.4 on 2024-07-11 12:31

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0014_alter_metadata_ent_date_alter_metadata_order_date_and_more'),
    ]

    operations = [
        migrations.CreateModel(
            name='CDSPositions',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('gen_bank_id', models.TextField()),
                ('cds_start', models.IntegerField()),
                ('cds_end', models.IntegerField()),
            ],
        ),
    ]
