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
import glob
import os
import re
import sys
import traceback
from datetime import datetime

from django.core.management.base import BaseCommand

import config.get_config as config
from core.helpers import read_csv_file, txt_to_list, compute_checksum
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
from helpers.color_log import logger

EMPTY_SAMPLES_GLOB = "*empty_samples*.txt"
META_DATA_GLOB = "*metadata*.csv"
SAMPLE_DIR_GLOB = "sample_*"
COUNT_TABLE_GLOB = "*count_table*.tsv"
PIPELINE_VIRSION_GLOB = "*pipeline_version*.txt"
CONTROL_DIR_GLOB = "sample_*-KO{neg,pos}"



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
                logger.warning(f"File type for {basename} not found")

            file, _ = File.objects.update_or_create(
                path=filepath,
                checksum=checksum,
                type=matching_file_type,
                sample=sample,
                plate=plate,
            )
            logger.info(f"Saved information about the file {file}")

        except FileNotFoundError as ex:
            logger.error(f"File with path {filepath} does not exist")
            logger.error(ex)

    def __process_sample_files(self, sample_dir, sample):
        file_types = FileType.objects.all()
        for file_type in file_types:
            file = glob.glob(os.path.join(sample_dir, f"*{file_type.postfix}*"))
            if file:
                self.__process_file(file[0], sample=sample)

    def __locate_plate_files(self, import_dir):
        print("----------------------------Import dir content--------------------------------------",
              os.listdir(import_dir))
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
                pseudonymized_id=item[columns.pseudonymized_id],
                plate=plate,
            )
            self.sample_id_dict[
                item[columns.pseudonymized_id]
            ] = False  # later on, we will check, if all samples have been imported
            self.__create_metadata(item, plate, well, sample)
        return plate

    """
    ------------------------------------- DEALING WITH COUNTS ----------------------------------------
    """

    def counts(self, plate, columns, counts_file, control=False, control_sample_id=None, control_type=None):
        """
        Process the counts file to create SampleCount object.
        :param counts_file: Path to the counts file
        :param columns: Columns to extract from the counts file
        :param plate: Plate object representing the processed plate
        :param control: Boolean indicating if the counts are for a control sample
        :param control_sample_is: The control sample id
        :param control_type: The control type
        :return: None
        """
        if control:
            Sample.objects.update_or_create(pseudonymized_id=control_sample_id,
                                            control=True,
                                            control_type=control_type)
            sample_id = control_sample_id
        else:
            sample_id = os.path.basename(counts_file).split("_")[0]
            if sample_id not in self.sample_id_dict:
                raise ValueError(f"Sample {sample_id} does not exist in the metadata file")
        data = read_csv_file(counts_file)
        try:
            sample = Sample.objects.get(pseudonymized_id=sample_id)

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
                        rpkm=float(item[columns.rpkm]) if item[columns.rpkm] != '' else None,
                        rpkm_proportions=float(item[columns.rpkm_proportions]) if item[columns.rpkm] != '' else None,
                        outlier=outlier,
                        coverage_threshold=float(item[columns.coverage_threshold]),
                        coverage=float(item[columns.coverage]),
                        coverage_status=item[columns.coverage_status],
                        readnum_status=item[columns.readnum_status],
                        readnum_threshold=float(item[columns.readnum_threshold]),
                        percentile_threshold=item[columns.percentile_threshold],
                        tax_id=int(item[columns.tax_id]),
                        scientific_name=item[columns.scientific_name],
                        DP=item[columns.DP] if columns.DP in item else None,
                        consensus_number_n=item[columns.consensus_number_n] if columns.consensus_number_n in item else None,
                        consensus_fraction_n=item[columns.consensus_fraction_n] if columns.consensus_fraction_n in item else None,
                    )
                    self.sample_id_dict[sample_id] = True
                    if _:
                        logger.info(
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
        empty_samples_message = "No empty samples found"
        if len(empty_samples) > 0:
            empty_samples_message = "Empty samples: " + ", ".join(empty_samples)
        logger.warning(empty_samples_message)
        if not all(self.sample_id_dict.values()):
            logger.error(
                "Not all samples have been imported. Please check the metadata file."
            )
            for key, value in self.sample_id_dict.items():
                if not value and key not in empty_samples:
                    logger.warning(f"Sample {key} has not been imported.")

    """
     ------------------------------------- MAIN -------------------------------------------------------------------
    """

    def handle(self, *args, **options):
        if len(sys.argv) < 2:
            logger.info("Usage: python import.py -d <import_dir>")
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
                all_dirs = glob.glob(os.path.join(import_dir, SAMPLE_DIR_GLOB))
                sample_dirs = []
                control_dirs = []
                for d in all_dirs:
                    if d.endswith('KOneg') or d.endswith('KOpos'):
                        control_dirs.append(d)
                    else:
                        sample_dirs.append(d)
                for sample_dir in sample_dirs:
                    logger.debug(
                        f"---------------- IMPORTING COUNTS FOR {sample_dir} ---------------- \n"
                    )
                    current_count_table_file = glob.glob(
                        os.path.join(sample_dir, COUNT_TABLE_GLOB)
                    )[0]
                    self.counts(plate, columns, current_count_table_file)
                    sample_id = os.path.basename(sample_dir).split("_")[1]
                    sample = Sample.objects.get(pseudonymized_id=sample_id)
                    logger.debug(
                        f"---------------- IMPORTING SAMPLE-RELATED FILES FOR THE SAMPLE {sample_dir} ---------------- \n"
                    )
                    self.__process_sample_files(sample_dir, sample)

                for control_dir in control_dirs:
                    logger.debug(f"----------------------- IMPORTING CONTROLS {control_dir} -----------------------\n")
                    if control_dir.endswith('pos'):
                        control_type = 'pos'
                        control_sample_id = f"control_pos_{plate.barcode}"
                    else:
                        control_sample_id = f"control_neg_{plate.barcode}"
                        control_type = 'neg'
                    current_count_table_file = glob.glob(
                        os.path.join(control_dir, COUNT_TABLE_GLOB)
                    )[0]
                    self.counts(plate, columns, current_count_table_file, control=True,
                                control_sample_id=control_sample_id, control_type=control_type)
                    sample = Sample.objects.get(pseudonymized_id=control_sample_id, control=True,
                                                control_type=control_type)


                    self.__process_sample_files(control_dir, sample)

                logger.debug(
                    f"\n --------------- IMPORTING PLATE-RELATED FILES FOR THE PLATE {plate} -------------------- \n"
                )
                self.__process_file(empty_samples_file, plate=plate)
                self.__process_file(metadata_file, plate=plate)
                self.__process_file(pipline_version_file, plate=plate)

                logger.debug(
                    "\n\n --------  CHECKING FOR MISSING SAMPLES...  ------------ \n"
                )
                self.check_imported_samples(empty_samples_file)

        except Exception as ex:
            logger.error(ex)
            traceback.print_exc()
