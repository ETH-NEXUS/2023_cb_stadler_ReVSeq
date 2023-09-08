"""
Import logic:

 ---Dealing with Files---
   - Verifies if given files exist and computes their checksum.
   - Categorizes and processes different files based on their postfix.
   - Locates and processes plate-specific files.

 ---Dealing with Metadata---
   - Processes metadata to extract panels.
   - Creates `Metadata`, `Plate`, `Well`, and `Sample` objects from given metadata.
   - Each metadata entry is associated with a plate, sample, and well.

 ---Dealing with Counts---
   - Processes count data to create `SampleCount` objects.
   - Each count entry is associated with a plate, sample, and substrain.

 ---Validation---
   - Checks if all samples have been successfully imported. If there are samples missing, the script identifies them.

 ---Main Execution---
   - Locates plate files in the provided import directory.
   - Processes metadata and count files.
   - Imports sample-related files and performs necessary validations.
   - If any error occurs, it logs the error and prints the traceback.
"""
import os
from datetime import datetime
import re
from django.core.management.base import BaseCommand
import traceback
from core.models import (
    SampleCount,
    Sample,
    Substrain,
    Plate,
    Well,
    Metadata,
    Panel,
    File,
    FileType,
)
import logging
import config.get_config as config
import sys
import glob
from core.helpers import read_csv_file, txt_to_list, compute_checksum

# from colorful_logger import logger as log

logger = logging.getLogger(__name__)

EMPTY_SAMPLES_GLOB = "*empty_samples*.txt"
META_DATA_GLOB = "*metadata*.csv"
SAMPLE_DIR_GLOB = "sample_*"
COUNT_TABLE_GLOB = "*count_table*.tsv"
PIPELINE_VIRSION_GLOB = "*pipeline_version*.txt"


class Command(BaseCommand):
    sample_id_dict = {}
    stats = {
        "processed_samples": 0,
        "errors": [],
        "warnings": [],
        "traces": [],
        "already_existed": 0,
        "overwritten": 0,
    }

    """
    ------------------------------------- ARGUMENTS ----------------------------------------
    """

    def add_arguments(self, parser):
        parser.add_argument(
            "--import_dir", "-d", type=str, help="Directory to import from"
        )

    """
    ------------------------------------- DEALING WITH FILES ----------------------------------------
    """

    def __process_file(self, filepath, sample=None, plate=None):
        try:
            if not os.path.exists(filepath):
                raise FileNotFoundError(f"File {filepath} does not exist")
            checksum = compute_checksum(filepath)
            basename = os.path.basename(filepath)
            file_types = FileType.objects.all()
            matching_file_type = None
            for ft in file_types:
                if ft.postfix in basename:
                    matching_file_type = ft
                    break
            if not matching_file_type:
                print(f"File type for {basename} not found")

            file, _ = File.objects.update_or_create(
                path=filepath,
                checksum=checksum,
                type=matching_file_type,
                sample=sample,
                plate=plate,
            )
            print(f"Saved information about the file {file}")

        except FileNotFoundError as ex:
            print(f"File with path {filepath} does not exist")
            logger.error(ex)

    def __process_sample_files(self, sample_dir, sample):
        file_types = FileType.objects.all()
        for file_type in file_types:
            file = glob.glob(os.path.join(sample_dir, f"*{file_type.postfix}*"))
            if file:
                self.__process_file(file[0], sample=sample)

    def __locate_plate_files(self, import_dir):
        metadata_file = glob.glob(os.path.join(import_dir, META_DATA_GLOB))[0]
        empty_samples_file = glob.glob(os.path.join(import_dir, EMPTY_SAMPLES_GLOB))[0]
        pipline_version_file = glob.glob(
            os.path.join(import_dir, PIPELINE_VIRSION_GLOB)
        )[0]
        if metadata_file and empty_samples_file:
            return metadata_file, empty_samples_file, pipline_version_file
        else:
            raise ValueError("Could not locate metadata and empty_samples files")

    """
    ------------------------------------- DEALING WITH METADATA ----------------------------------------
    """

    def __extract_panels(self, item):
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

    def __create_metadata(self, item, plate, well, sample):
        order_date = datetime.strptime(item["Order date"], "%Y-%m-%d")
        ent_date = datetime.strptime(item["ent date"], "%Y-%m-%d")
        treatment_type = item["treatment_type"] if item["treatment_type"] else None
        data = self.__extract_panels(item)
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
            self.__create_metadata(item, plate, well, sample)
        return plate

    """
    ------------------------------------- DEALING WITH COUNTS ----------------------------------------
    """

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

    """
    ------------------------------------- CHECKING IF ALL SAMPLES HAVE BEEN IMPORTED ----------------------------------------
    """

    def check_imported_samples(self, file_name):
        empty_samples = txt_to_list(file_name)
        print("Empty samples:", ", ".join(empty_samples).replace("\n", ""))
        if not all(self.sample_id_dict.values()):
            print("Not all samples have been imported. Please check the metadata file.")
            for key, value in self.sample_id_dict.items():
                if not value and key not in empty_samples:
                    print(f"Sample {key} has not been imported.")

    """
     ------------------------------------- MAIN -------------------------------------------------------------------
    """

    def handle(self, *args, **options):
        if len(sys.argv) < 2:
            print("Usage: python import.py -d <import_dir>")
            sys.exit(1)
        import_dir = options.get("import_dir")
        (
            metadata_file,
            empty_samples_file,
            pipline_version_file,
        ) = self.__locate_plate_files(import_dir)

        try:
            if metadata_file and empty_samples_file:
                columns = config.read_config("import_columns")
                plate = self.metadata(columns, metadata_file)

                sample_dirs = glob.glob(os.path.join(import_dir, SAMPLE_DIR_GLOB))
                for sample_dir in sample_dirs:
                    print(
                        f"---------------- IMPORTING COUNTS FOR {sample_dir} ---------------- \n"
                    )
                    current_count_table_file = glob.glob(
                        os.path.join(sample_dir, COUNT_TABLE_GLOB)
                    )[0]
                    self.counts(plate, columns, current_count_table_file)
                    sample_id = os.path.basename(sample_dir).split("_")[1]
                    sample = Sample.objects.get(pseudoanonymized_id=sample_id)
                    print(
                        f"---------------- IMPORTING SAMPLE-RELATED FILES FOR THE SAMPLE {sample_dir} ---------------- \n"
                    )
                    self.__process_sample_files(sample_dir, sample)

                print(
                    f"\n --------------- IMPORTING PLATE-RELATED FILES FOR THE PLATE {plate} -------------------- \n"
                )
                self.__process_file(empty_samples_file, plate=plate)
                self.__process_file(metadata_file, plate=plate)
                self.__process_file(pipline_version_file, plate=plate)

                print("\n\n --------  CHECKING FOR MISSING SAMPLES...  ------------ \n")
                self.check_imported_samples(empty_samples_file)

        except Exception as ex:
            logger.error(ex)
            traceback.print_exc()
