# Generated by Django 4.2.4 on 2024-02-05 09:59

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0005_alter_samplecount_not_valid'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='samplecount',
            name='not_valid',
        ),
        migrations.AddField(
            model_name='sample',
            name='not_valid',
            field=models.BooleanField(blank=True, default=False, null=True),
        ),
    ]
