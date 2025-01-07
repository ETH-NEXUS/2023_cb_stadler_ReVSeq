from django.core.management.base import BaseCommand
from core.models import Sample, Plate
from helpers.color_log import logger

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('mode', type=str, choices=['metagenomics', 'alignment'], help='Mode to set for samples')
        parser.add_argument('--add_prefix', action='store_true', help='Whether to add a prefix to the pseudonymized_id')
        parser.add_argument('--prefix', type=str, default='', help='Prefix to add to the pseudonymized_id if add_prefix is set')
        parser.add_argument('--plate', type=str, help='Plate of the affected plate')

    def handle(self, *args, **options):
        mode = options['mode']
        add_prefix = options['add_prefix']
        prefix = options['prefix']
        plate_barcode = options['plate']

        plate = Plate.objects.get(barcode=plate_barcode)
        if plate:
            samples = Sample.objects.filter(plate=plate)
            for sample in samples:
                sample.mode = mode
                if add_prefix:
                    sample.pseudonymized_id = prefix + '-' +  sample.pseudonymized_id
                sample.save()
                logger.info(f"Mode for {sample} changed to {mode}")
        else:
            logger.error(f"Plate with barcode {plate_barcode} not found")

# ./manage.py change_mode alignment --add_prefix --prefix a' --plate <PlateBarcode>