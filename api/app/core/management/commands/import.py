"""
Import logic:
1) First we create a new plate object based on the barcode in the name of the metadata file
2) We create a metadata object for every well (every row represents a well) in the metadata file. 
Also, for every well we create the well object and the sample object.
3)The sample counts file refers to the whole plate. We create it at the end. 

"""
import os
from datetime import datetime
import re
from django.core.management.base import BaseCommand
import traceback
import csv
from core.models import SampleCount, Sample, Substrain, Plate, Well, Metadata, Panel
import logging
import config.get_config as config

logger = logging.getLogger(__name__)


def read_csv_file(input_file):
    with open(input_file, "r", encoding="utf-8-sig") as f:
        dialect = csv.Sniffer().sniff(f.readline())
        f.seek(0)
        return list(csv.DictReader(f, dialect=dialect))


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--metadata_file", "-m", type=str, help="File with the metadata to upload"
        )
        parser.add_argument(
            "--counts_file", "-c", type=str, help="File with counts to upload"
        )

    def extract_panels(self, item):
        data = []
        for key, value in item.items():
            match = re.match(r"([a-zA-Z0-9]+) \(\d+\)", key)
            if match:
                relevant_key = match.group(1).strip()
                panel = Panel.objects.get(name=relevant_key)
                strain = panel.strain
                data.append(
                    {
                        "panel": match.group(1),
                        "value": value,
                        "strain": strain.name,
                    }
                )
        return data

    def create_metadata(self, item, plate, well, sample):
        order_date = datetime.strptime(item["Order date"], "%Y-%m-%d")
        ent_date = datetime.strptime(item["ent date"], "%Y-%m-%d")
        treatment_type = item["treatment_type"] if item["treatment_type"] else None
        data = self.extract_panels(item)
        return Metadata.objects.create(
            plate=plate,
            sample=sample,
            well=well,
            prescriber=item["Prescriber canton"],
            order_date=order_date,
            ent_date=ent_date,
            treatment_type=treatment_type,
            data=data,
        )

    def metadata(self, columns, **options):
        """
        Process the metadata file to create Plate, Well, Sample, and Metadata objects.
        :param options: Command-line options including metadata file path
        :return: Plate object representing the processed plate
        """

        barcode = os.path.basename(options.get("metadata_file")).split("_")[0]
        plate, _ = Plate.objects.get_or_create(barcode=barcode)
        data = read_csv_file(options.get("metadata_file"))

        for item in data:
            well, _ = Well.objects.get_or_create(
                plate=plate, location=item[columns.deep_well_location]
            )
            sample, _ = Sample.objects.get_or_create(
                well=well,
                sample_number=item[columns.sample_number],
                pseudoanonymized_id=item[columns.pseudoanonymized_id],
                plate=plate,
            )
            self.create_metadata(item, plate, well, sample)
        return plate

    def counts(self, plate, columns, **options):
        """
        Process the counts file to create SampleCount object.
        :param options: Command-line options including counts file path
        :return: None
        """
        file_name = options.get("counts_file")
        data = read_csv_file(file_name)
        sample_id = os.path.basename(file_name).split("_")[0]
        try:
            sample = Sample.objects.get(pseudoanonymized_id=sample_id)
            for item in data:
                substrain_name = item[columns.name]
                outlier = False
                if item[columns.outlier] == "*":
                    outlier = True
                try:
                    substrain = Substrain.objects.get(name=substrain_name)
                    sampleCount, _ = SampleCount.objects.update_or_create(
                        plate=plate,
                        sample=sample,
                        substrain=substrain,
                        aligned=int(item[columns.aligned]),
                        length=int(item[columns.length]),
                        rpkm=float(item[columns.rpkm]),
                        rpkm_proportions=float(item[columns.rpkm_proportions]),
                        normcounts=float(item[columns.normcounts]),
                        outlier=outlier,
                    )
                except Substrain.DoesNotExist:
                    logger.error(f"Substrain {substrain_name} does not exist")
                    raise ValueError(f"Substrain {substrain_name} does not exist")
        except Sample.DoesNotExist:
            logger.error(f"Sample {sample_id} does not exist")
            raise ValueError(f"Sample {sample_id} does not exist")

    def handle(self, *args, **options):
        try:
            columns = config.read_config("import_columns")
            plate = self.metadata(columns, **options)
            self.counts(plate, columns, **options)
            print("Import successful")
        except Exception as ex:
            logger.error(ex)
            traceback.print_exc()
