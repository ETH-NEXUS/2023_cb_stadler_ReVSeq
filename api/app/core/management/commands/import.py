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
import sys
import glob

# from colorful_logger import logger as log

logger = logging.getLogger(__name__)


def list_files(startpath):
    for root, dirs, files in os.walk(startpath):
        level = root.replace(startpath, "").count(os.sep)
        indent = "│   " * (level - 1) if level > 0 else ""
        dir_basename = os.path.basename(root)
        if root != startpath:
            print("{}├── {}/".format(indent, dir_basename))
        subindent = "│   " * level
        file_indent = "{}├── ".format(subindent)
        for i, f in enumerate(files):
            if i == len(files) - 1 and not dirs:
                file_indent = "{}└── ".format(subindent)
            print("{}{}".format(file_indent, f))


def txt_to_list(input_file):
    with open(input_file, "r", encoding="utf-8-sig") as f:
        return f.readlines()


def read_csv_file(input_file):
    with open(input_file, "r", encoding="utf-8-sig") as f:
        dialect = csv.Sniffer().sniff(f.readline())
        f.seek(0)
        return list(csv.DictReader(f, dialect=dialect))


class Command(BaseCommand):
    sample_id_dict = {}

    def add_arguments(self, parser):
        parser.add_argument(
            "--import_dir", "-d", type=str, help="Directory to import from"
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
        metadata, _ = Metadata.objects.update_or_create(
            plate=plate,
            sample=sample,
            well=well,
            prescriber=item["Prescriber canton"],
            order_date=order_date,
            ent_date=ent_date,
            treatment_type=treatment_type,
            data=data,
        )
        return metadata

    def metadata(self, columns, metadata_file):
        """
        Process the metadata file to create Plate, Well, Sample, and Metadata objects.
        :param metadata_file: Path to the metadata file
        :return: Plate object representing the processed plate
        """

        barcode = os.path.basename(metadata_file.split("_")[0])
        plate, _ = Plate.objects.get_or_create(barcode=barcode)
        data = read_csv_file(metadata_file)

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
            self.sample_id_dict[
                item[columns.pseudoanonymized_id]
            ] = False  # later on, we will check, if all samples have been imported
            self.create_metadata(item, plate, well, sample)
        return plate

    def counts(self, plate, columns, counts_file):
        """
        Process the counts file to create SampleCount object.
        :param counts_file: Path to the counts file
        :return: None
        """
        data = read_csv_file(counts_file)
        sample_id = os.path.basename(counts_file).split("_")[0]
        if sample_id not in self.sample_id_dict:
            raise ValueError(f"Sample {sample_id} does not exist in the metadata file")
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
                        # normcounts=float(item[columns.normcounts]),
                        outlier=outlier,
                        qc_status=item[columns.qc_status],
                        coverage_threshold=float(item[columns.coverage_threshold]),
                        coverage=float(item[columns.coverage]),
                    )
                    self.sample_id_dict[sample_id] = True
                    if _:
                        print(
                            f"Imported sample counts for {substrain} in the sample {sample}"
                        )
                except Substrain.DoesNotExist:
                    logger.error(f"Substrain {substrain_name} does not exist")
                    raise ValueError(f"Substrain {substrain_name} does not exist")

        except Sample.DoesNotExist:
            logger.error(f"Sample {sample_id} does not exist")
            raise ValueError(f"Sample {sample_id} does not exist")

    def check_imported_samples(self, file_name):
        empty_samples = txt_to_list(file_name)
        print("Empty samples:", ", ".join(empty_samples).replace("\n", ""))
        if not all(self.sample_id_dict.values()):
            print("Not all samples have been imported. Please check the metadata file.")
            for key, value in self.sample_id_dict.items():
                if not value and key not in empty_samples:
                    print(f"Sample {key} has not been imported.")
            raise ValueError(
                "Not all samples have been imported. Please check the metadata file."
            )

    def __locate_files(self, import_dir):
        metadata_file = glob.glob(os.path.join(import_dir, "*metadata*.csv"))[0]
        empty_samples_file = glob.glob(os.path.join(import_dir, "*empty_samples*.txt"))[
            0
        ]
        if metadata_file and empty_samples_file:
            return metadata_file, empty_samples_file
        else:
            raise ValueError("Could not locate metadata and empty_samples files")

    def handle(self, *args, **options):
        if len(sys.argv) < 2:
            print("Usage: python import.py -d <import_dir>")
            sys.exit(1)
        import_dir = options.get("import_dir")
        metadata_file, empty_samples_file = self.__locate_files(import_dir)

        try:
            if metadata_file and empty_samples_file:
                columns = config.read_config("import_columns")
                plate = self.metadata(columns, metadata_file)
                count_files = glob.glob(
                    os.path.join(import_dir, "sample_*", "*count_table*.tsv")
                )
                for file in count_files:
                    self.counts(plate, columns, file)
                    print(f"Imported {file}")
                self.check_imported_samples(empty_samples_file)
                print("Import successful")

        except Exception as ex:
            logger.error(ex)
            traceback.print_exc()
