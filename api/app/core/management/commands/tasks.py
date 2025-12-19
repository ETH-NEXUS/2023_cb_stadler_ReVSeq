from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Optional

from django.core.management.base import BaseCommand
from django.db import transaction

from helpers.color_log import logger
from core.models import Sample, File, FileType


# ---------------------------------------------------------------------
# Hardcoded input data (pseudonymized_id, folder_path)
# ---------------------------------------------------------------------
SAMPLES_AND_FOLDERS: list[tuple[str, str]] = [
    ("m2-Dq37e8", "/data/RVSeqPlate10-m2/sample_m2-Dq37e8/"),
    ("m2-wfEdma", "/data/RVSeqPlate10-m2/sample_m2-wfEdma/"),
    ("m2-9eyPAv", "/data/RVSeqPlate10-m2/sample_m2-9eyPAv/"),
    ("m2-GRjNzo", "/data/RVSeqPlate11-m2/sample_m2-GRjNzo/"),
    ("m2-Kg3Fdm", "/data/RVSeqPlate11-m2/sample_m2-Kg3Fdm/"),
    ("m2-SrqwxQ", "/data/RVSeqPlate11-m2/sample_m2-SrqwxQ/"),
    ("m2-mzMRtV", "/data/RVSeqPlate13-m2/sample_m2-mzMRtV/"),
    ("m2-zwqknf", "/data/RVSeqPlate3-m2/sample_m2-zwqknf/"),
    ("m2-3oZJ5x", "/data/RVSeqPlate3-m2/sample_m2-3oZJ5x/"),
    ("m2-43Vmvj", "/data/RVSeqPlate3-m2/sample_m2-43Vmvj/"),
    ("m2-NmHGDa", "/data/RVSeqPlate4-m2/sample_m2-NmHGDa/"),
    ("m2-u89dBs", "/data/RVSeqPlate4-m2/sample_m2-u89dBs/"),
    ("m2-gnDemF", "/data/RVSeqPlate4-m2/sample_m2-gnDemF/"),
    ("m2-HXHQAT", "/data/RVSeqPlate4-m2/sample_m2-HXHQAT/"),
    ("m2-RVSeqPlate5-KOpos", "/data/RVSeqPlate5-m2/sample_m2-RVSeqPlate5-KOpos/"),
    ("m2-w9v3Qv", "/data/RVSeqPlate5-m2/sample_m2-w9v3Qv/"),
    ("m2-F2HifQ", "/data/RVSeqPlate5-m2/sample_m2-F2HifQ/"),
    ("m2-EJ7fEz", "/data/RVSeqPlate5-m2/sample_m2-EJ7fEz/"),
    ("m2-RVSeqPlate6-KOpos", "/data/RVSeqPlate6-m2/sample_m2-RVSeqPlate6-KOpos/"),
    ("m2-4L3pXZ", "/data/RVSeqPlate6-m2/sample_m2-4L3pXZ/"),
    ("m2-NDEjHn", "/data/RVSeqPlate7-m2/sample_m2-NDEjHn/"),
    ("m2-Y4rZ7q", "/data/RVSeqPlate8-m2/sample_m2-Y4rZ7q/"),
    ("m2-XHdiUs", "/data/RVSeqPlate8-m2/sample_m2-XHdiUs/"),
]


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def md5sum(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def detect_postfix(filename: str) -> str:
    """
    FileType.postfix should match your pipeline expectations.
    Path.suffix would only yield ".gz", so we detect multi-suffix endings.
    """
    lowered = filename.lower()
    for ending in (".txt.gz",):
        if lowered.endswith(ending):
            return ending
    return Path(filename).suffix or lowered


def pick_major_minor_chr_files(folder: Path) -> tuple[Optional[Path], Optional[Path]]:
    """
    Pick chr_file_major.txt.gz and chr_file_minor.txt.gz if present.
    """
    major_chr = folder / "chr_file_major.txt.gz"
    minor_chr = folder / "chr_file_minor.txt.gz"

    major = major_chr if major_chr.exists() else None
    minor = minor_chr if minor_chr.exists() else None
    return major, minor


def upsert_file(*, sample: Sample, file_path: Path, dry_run: bool) -> None:
    postfix = detect_postfix(file_path.name)

    if dry_run:
        logger.info(f"[DRY RUN] Would attach file {file_path} to sample {sample} (postfix={postfix})")
        return

    checksum = md5sum(file_path)
    ft, _ = FileType.objects.get_or_create(postfix=postfix)

    obj, created = File.objects.get_or_create(
        path=str(file_path),
        defaults={
            "checksum": checksum,
            "type": ft,
            "sample": sample,
            "plate": sample.plate,
        },
    )

    if created:
        logger.info(
            f"Attached NEW File {file_path} to sample {sample} "
            f"(checksum={checksum}, postfix={postfix})"
        )
        return

    changed = False
    if obj.checksum != checksum:
        obj.checksum = checksum
        changed = True
    if obj.type_id != ft.id:
        obj.type = ft
        changed = True
    if obj.sample_id != sample.id:
        obj.sample = sample
        changed = True
    if obj.plate_id != sample.plate_id:
        obj.plate = sample.plate
        changed = True

    if changed:
        obj.save()
        logger.info(
            f"Updated existing File {file_path} for sample {sample} "
            f"(checksum={checksum}, postfix={postfix})"
        )
    else:
        logger.debug(f"No changes for existing File {file_path} (already attached to sample {sample})")


# ---------------------------------------------------------------------
# Command
# ---------------------------------------------------------------------
class Command(BaseCommand):
    help = "Attach chr_file_major.txt.gz / chr_file_minor.txt.gz to Samples via File model."

    def add_arguments(self, parser):
        parser.add_argument(
            "--dry-run",
            action="store_true",
            default=False,
            help="Do not write to DB, only log what would be attached.",
        )

    def handle(self, *args, **options):
        dry_run: bool = options["dry_run"]

        logger.info(f"Starting attach_chr_files (dry_run={dry_run})")
        logger.info(f"Processing {len(SAMPLES_AND_FOLDERS)} sample folder entries")

        for pseudonymized_id, folder_str in SAMPLES_AND_FOLDERS:
            folder = Path(folder_str)

            sample = Sample.objects.filter(pseudonymized_id=pseudonymized_id).select_related("plate").first()
            if sample is None:
                logger.warning(f"Sample not found in DB: {pseudonymized_id} (folder={folder})")
                continue

            if not folder.exists() or not folder.is_dir():
                logger.warning(f"Folder missing or not a directory for sample {sample}: {folder}")
                continue

            major_chr, minor_chr = pick_major_minor_chr_files(folder)

            if major_chr is None and minor_chr is None:
                logger.warning(
                    f"Sample {sample}: neither chr_file_major.txt.gz nor chr_file_minor.txt.gz found in {folder}"
                )
                continue

            if major_chr is None:
                logger.warning(f"Sample {sample}: chr_file_major.txt.gz missing in {folder}; will attach only MINOR if present")
            if minor_chr is None:
                logger.warning(f"Sample {sample}: chr_file_minor.txt.gz missing in {folder}; will attach only MAJOR if present")

            with transaction.atomic():
                if major_chr is not None:
                    upsert_file(sample=sample, file_path=major_chr, dry_run=dry_run)
                if minor_chr is not None:
                    upsert_file(sample=sample, file_path=minor_chr, dry_run=dry_run)

        logger.info("Done.")