import csv
from django.core.management.base import BaseCommand
from django.db import transaction
from core.models import Sample, File, FileType
from core.helpers import compute_checksum
from helpers.color_log import logger
import os

class Command(BaseCommand):
    help = "Import consensus .gz files and .bed files from a CSV"
    def add_arguments(self, parser):
        parser.add_argument(
            "csv_path",
            type=str,
            help="Path to consensus_files.csv",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Do not save anything, only print actions",
        )

    def handle(self, *args, **options):
        csv_path = options["csv_path"]
        dry_run = options["dry_run"]
        logger.info(f"Reading CSV: {csv_path}")
        created = 0
        skipped = 0
        errors = 0

        with transaction.atomic():
            with open(csv_path, newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    pseudonymized_id = row["id"].strip()
                    gz_path = row["path"].strip()
                    bed_path = row["bed"].strip()

                    try:
                        sample = Sample.objects.get(pseudonymized_id=pseudonymized_id)
                    except Sample.DoesNotExist:
                        logger.error(f"Sample {pseudonymized_id} not found.")
                        errors += 1
                        continue
                    plate = sample.plate
                    # Process both files
                    for file_path in (gz_path, bed_path):
                        postfix = self._detect_filetype(file_path)

                        if postfix is None:
                            logger.error(
                                f"Cannot determine file type for {file_path}"
                            )
                            errors += 1
                            continue

                        file_type, _ = FileType.objects.get_or_create(postfix=postfix)

                        if File.objects.filter(path=file_path).exists():
                            logger.info(f"SKIP: {file_path} already exists")
                            skipped += 1
                            continue

                        logger.info(f"ADD: {file_path} (type={postfix}) "
                                    f"→ sample={sample.pseudonymized_id}, plate={plate}")
                        if dry_run:
                            continue
                        checksum = compute_checksum(file_path)
                        File.objects.create(
                            path=file_path,
                            checksum=checksum,
                            type=file_type,
                            sample=sample,
                            plate=plate,
                        )

                        created += 1

            if dry_run:
                logger.info("Dry run enabled → rolling back transaction")
                raise transaction.TransactionManagementError("Dry run: rollback")

        logger.info("=== Import summary ===")
        logger.info(f"Created: {created}")
        logger.info(f"Skipped: {skipped}")
        logger.info(f"Errors: {errors}")

    def _detect_filetype(self, path: str):
        """
        Match path suffix to FileType.postfix.
        Handles things like:
        - .fa
        - .fa.gz
        - .txt
        - .txt.gz
        - .gz
        - .bed
        """

        filename = os.path.basename(path)

        # Try longest postfixes first
        candidate_postfixes = FileType.objects.values_list("postfix", flat=True)
        candidate_postfixes = sorted(candidate_postfixes, key=len, reverse=True)

        for postfix in candidate_postfixes:
            if filename.endswith(postfix):
                return postfix

        return None