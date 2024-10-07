from django.core.management.base import BaseCommand
from core.models import Sample


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('mode', type=str, choices=['metagenomics', 'alignment'], help='Mode to set for samples')
        parser.add_argument('--add_prefix', action='store_true', help='Whether to add a prefix to the pseudonymized_id')
        parser.add_argument('--prefix', type=str, default='', help='Prefix to add to the pseudonymized_id if add_prefix is set')

    def handle(self, *args, **options):
        mode = options['mode']
        add_prefix = options['add_prefix']
        prefix = options['prefix']

        samples = Sample.objects.all()
        for sample in samples:
            sample.mode = mode
            if add_prefix:
                sample.pseudonymized_id = prefix + '-' +  sample.pseudonymized_id
            sample.save()
            print(f"Mode for {sample} changed to {mode}")