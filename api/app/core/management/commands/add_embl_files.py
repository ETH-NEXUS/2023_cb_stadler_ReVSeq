from django.core.management.base import BaseCommand
from core.models import Sample, File, Plate, FileType
from helpers.color_log import logger
from core.helpers import compute_checksum

file_data = [
    {"plate": "RVSeqPlate4-m2", "sample": "RyXauM", "file": "/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-RyXauM/RyXauM-m2.embl.gz"},

]

class Command(BaseCommand):
    help = 'Imports file paths and assigns them to samples from a text file.'

    def handle(self, *args, **kwargs):
        for entry in file_data:
            plate_barcode = entry["plate"]
            sample_id = entry["sample"]
            file_path = entry["file"]

            try:
                plate = Plate.objects.get(barcode=plate_barcode)
                sample = Sample.objects.get(pseudonymized_id=sample_id)
                file_type, created = FileType.objects.get_or_create(postfix="embl")
                checksum = compute_checksum(file_path)

                file_instance = File.objects.create(
                    path=file_path,
                    checksum=checksum,
                    type=file_type,
                    sample=sample,
                    plate=plate
                )
                logger.info(f'File {file_path} created with type {file_type.postfix} for sample {sample_id} on plate {plate_barcode}')
                logger.info(f'Added file {file_path} for sample {sample_id} on plate {plate_barcode}')
            except (Plate.DoesNotExist, Sample.DoesNotExist) as e:
                logger.error(f'Error adding file {file_path} for sample {sample_id} on plate {plate_barcode}: {e}')