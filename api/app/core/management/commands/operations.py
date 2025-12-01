import os
import shutil
import hashlib

from django.core.management.base import BaseCommand
from django.db import IntegrityError, transaction

from core.models import Sample, File, FileType
from helpers.color_log import logger  # if you prefer, you can also use self.stdout.write


CONSENSUS_SUFFIX = "consensus_upload.gz"
NEW_SUFFIX = "consensus_upload.fasta.gz"
NEW_FILETYPE_POSTFIX = ".fasta.gz"


def compute_md5(path: str, chunk_size: int = 8192) -> str:
    """
    Compute MD5 checksum of a file at the given path.
    """
    md5 = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            md5.update(chunk)
    return md5.hexdigest()


class Command(BaseCommand):
    help = (
        "For each Sample, find File entries whose path ends with 'consensus_upload.gz', "
        "copy the file to the same directory with suffix 'consensus_upload.fasta.gz', "
        "and create a new File entry with the new path and checksum."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "--dry-run",
            action="store_true",
            dest="dry_run",
            help="Do not modify files or database, just log what would be done.",
        )

    def handle(self, *args, **options):
        dry_run = options.get("dry_run", False)

        # Ensure we have a FileType for the new suffix
        if dry_run:
            logger.info(
                "DRY RUN: would ensure FileType with postfix '%s' exists.",
                NEW_FILETYPE_POSTFIX,
            )
            filetype = None
        else:
            filetype, _ = FileType.objects.get_or_create(postfix=NEW_FILETYPE_POSTFIX)
            logger.info("Using FileType with postfix '%s' (id=%s).", filetype.postfix, filetype.id)

        samples = Sample.objects.all()
        logger.info("Processing %d samples.", samples.count())

        total_original_files = 0
        total_new_files_created = 0

        for sample in samples:
            # Find all File entries for this sample with the old consensus suffix
            files = File.objects.filter(sample=sample, path__endswith=CONSENSUS_SUFFIX)
            if not files.exists():
                continue

            logger.info(
                "Sample %s: found %d consensus files to process.",
                sample.pseudonymized_id,
                files.count(),
            )
            total_original_files += files.count()

            for original_file in files:
                old_path = original_file.path
                dir_name, base_name = os.path.split(old_path)

                if not base_name.endswith(CONSENSUS_SUFFIX):
                    # Defensive: should not happen due to the filter, but we check anyway.
                    logger.warning(
                        "File path %s does not actually end with %s; skipping.",
                        old_path,
                        CONSENSUS_SUFFIX,
                    )
                    continue

                new_base_name = base_name[: -len(CONSENSUS_SUFFIX)] + NEW_SUFFIX
                new_path = os.path.join(dir_name, new_base_name)

                logger.info(
                    "Sample %s: copying %s -> %s",
                    sample.pseudonymized_id,
                    old_path,
                    new_path,
                )

                # Check that the source exists
                if not os.path.exists(old_path):
                    logger.error("Source file does not exist on disk: %s", old_path)
                    continue

                # If the target already exists, we can either skip or just reuse it.
                if os.path.exists(new_path):
                    logger.warning("Target file already exists: %s", new_path)
                else:
                    if dry_run:
                        logger.info("DRY RUN: would copy file on disk.")
                    else:
                        try:
                            os.makedirs(dir_name, exist_ok=True)
                            shutil.copy2(old_path, new_path)
                        except Exception as e:
                            logger.error(
                                "Failed to copy '%s' to '%s': %s",
                                old_path,
                                new_path,
                                e,
                            )
                            continue

                if dry_run:
                    logger.info(
                        "DRY RUN: would compute checksum and create File entry for %s",
                        new_path,
                    )
                    continue

                # Compute checksum for the new file
                try:
                    checksum = compute_md5(new_path)
                except Exception as e:
                    logger.error("Failed to compute MD5 for %s: %s", new_path, e)
                    continue

                # Create the new File row
                try:
                    with transaction.atomic():
                        new_file = File.objects.create(
                            path=new_path,
                            checksum=checksum,
                            type=filetype,
                            sample=sample,
                            plate=sample.plate,
                        )
                    total_new_files_created += 1
                    logger.info(
                        "Created new File entry id=%s for path %s",
                        new_file.id,
                        new_path,
                    )
                except IntegrityError:
                    # path is unique, so this can happen if the DB already has this File
                    logger.warning(
                        "File entry with path %s already exists in DB; skipping create.",
                        new_path,
                    )
                except Exception as e:
                    logger.error(
                        "Unexpected error while creating File entry for %s: %s",
                        new_path,
                        e,
                    )

        logger.info(
            "Done. Processed %d original consensus files; created %d new .fasta.gz File entries%s.",
            total_original_files,
            total_new_files_created,
            " (dry run)" if dry_run else "",
        )



# from pathlib import Path
#
# from django.core.management.base import BaseCommand, CommandError
#
# from core.models import Sample
# from helpers.color_log import logger
#
#
# class Command(BaseCommand):
#     help = (
#         "Reset job_id and analysis_job_id to None for samples listed in a text file.\n"
#         "The file should contain one pseudonymized_id per line."
#     )
#
#     def add_arguments(self, parser):
#         parser.add_argument(
#             "file_path",
#             type=str,
#             help="Path to the text file with pseudonymized_ids (one per line).",
#         )
#
#     def handle(self, *args, **options):
#         file_path = options["file_path"]
#         path = Path(file_path)
#
#         if not path.exists():
#             raise CommandError(f"File not found: {path}")
#
#         # Read IDs, strip whitespace, ignore empty lines
#         ids = [
#             line.strip()
#             for line in path.read_text().splitlines()
#             if line.strip()
#         ]
#
#         if not ids:
#             self.stdout.write(self.style.WARNING("No IDs found in file. Nothing to do."))
#             return
#
#         self.stdout.write(f"Read {len(ids)} pseudonymized IDs from {path}.")
#
#         updated = 0
#         not_found = []
#
#         for pid in ids:
#             qs = Sample.objects.filter(pseudonymized_id=pid)
#             if not qs.exists():
#                 not_found.append(pid)
#                 logger.warning(f"No Sample found for pseudonymized_id={pid}")
#                 continue
#
#             count = qs.update(job_id=None, analysis_job_id=None)
#             updated += count
#             logger.info(f"Reset job_id and analysis_job_id for {count} sample(s) with pseudonymized_id={pid}")
#
#         self.stdout.write(self.style.SUCCESS(f"Updated {updated} Sample objects."))
#
#         if not_found:
#             self.stdout.write(
#                 self.style.WARNING(
#                     f"{len(not_found)} IDs had no matching Sample:\n" +
#                     "\n".join(not_found)
#                 )
#             )