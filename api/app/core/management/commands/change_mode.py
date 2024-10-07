from django.core.management.base import BaseCommand
from os import environ
from core.models import Sample


class Command(BaseCommand):

    def handle(self, *args, **options):
        samples = Sample.objects.all()
        for sample in samples:
            sample.mode = "alignment"
            sample.pseudonymized_id = "al-" + sample.pseudonymized_id
            sample.save()
            print(f"Mode for {sample} changed to alignment")