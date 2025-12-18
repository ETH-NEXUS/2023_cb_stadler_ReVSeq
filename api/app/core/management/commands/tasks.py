from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Optional, Tuple, List

from django.core.management.base import BaseCommand
from django.db import transaction

from helpers.color_log import logger
from core.models import Sample, File, FileType


# ---------------------------------------------------------------------
# Hardcoded input data (pseudonymized_id, folder_path)
# ---------------------------------------------------------------------
SAMPLES_AND_FOLDERS: list[tuple[str, str]] = [
    ("m2-Dq37e8", "/data/revseq/results/gather_results/RVSeqPlate10-m2/sample_m2-Dq37e8/"),
    ("m2-wfEdma", "/data/revseq/results/gather_results/RVSeqPlate10-m2/sample_m2-wfEdma/"),
    ("m2-9eyPAv", "/data/revseq/results/gather_results/RVSeqPlate10-m2/sample_m2-9eyPAv/"),
    ("m2-GRjNzo", "/data/revseq/results/gather_results/RVSeqPlate11-m2/sample_m2-GRjNzo/"),
    ("m2-Kg3Fdm", "/data/revseq/results/gather_results/RVSeqPlate11-m2/sample_m2-Kg3Fdm/"),
    ("m2-SrqwxQ", "/data/revseq/results/gather_results/RVSeqPlate11-m2/sample_m2-SrqwxQ/"),
    ("m2-mzMRtV", "/data/revseq/results/gather_results/RVSeqPlate13-m2/sample_m2-mzMRtV/"),
    ("m2-zwqknf", "/data/revseq/results/gather_results/RVSeqPlate3-m2/sample_m2-zwqknf/"),
    ("m2-3oZJ5x", "/data/revseq/results/gather_results/RVSeqPlate3-m2/sample_m2-3oZJ5x/"),
    ("m2-43Vmvj", "/data/revseq/results/gather_results/RVSeqPlate3-m2/sample_m2-43Vmvj/"),
    ("m2-NmHGDa", "/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-NmHGDa/"),
    ("m2-u89dBs", "/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-u89dBs/"),
    ("m2-gnDemF", "/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-gnDemF/"),
    ("m2-HXHQAT", "/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-HXHQAT/"),
    ("m2-RVSeqPlate5-KOpos", "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-RVSeqPlate5-KOpos/"),
    ("m2-w9v3Qv", "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-w9v3Qv/"),
    ("m2-F2HifQ", "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-F2HifQ/"),
    ("m2-EJ7fEz", "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-EJ7fEz/"),
    ("m2-RVSeqPlate6-KOpos", "/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-RVSeqPlate6-KOpos/"),
    ("m2-4L3pXZ", "/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-4L3pXZ/"),
    ("m2-NDEjHn", "/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-NDEjHn/"),
    ("m2-Y4rZ7q", "/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-Y4rZ7q/"),
    ("m2-XHdiUs", "/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-XHdiUs/"),
    # note: zwqknf appears twice in your list; we keep first and silently dedupe by File.path uniqueness
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
    for ending in (".embl.gz", ".fa.gz", ".fasta.gz"):
        if lowered.endswith(ending):
            return ending
    # fallback (keeps DB consistent even for unexpected names)
    return Path(filename).suffix or lowered


def pick_major_minor_files(folder: Path) -> tuple[Optional[Path], Optional[Path]]:
    """
    New rule: major/minor are independent:
      - prefer EMBL if present, else FASTA
    """
    # Major candidates
    major_embl = next(folder.glob("*_consensus_major.embl.gz"), None)
    major_fa = next(folder.glob("*_consensus_major.fa.gz"), None)

    # Minor candidates
    minor_embl = next(folder.glob("*_consensus_minor.embl.gz"), None)
    minor_fa = next(folder.glob("*_consensus_minor.fa.gz"), None)

    major = major_embl or major_fa
    minor = minor_embl or minor_fa
    return major, minor


def upsert_file(*, sample: Sample, file_path: Path, dry_run: bool) -> None:
    postfix = detect_postfix(file_path.name)

    if dry_run:
        logger.info(f"[DRY RUN] Would attach file {file_path} to sample {sample} (postfix={postfix})")
        return

    checksum = md5sum(file_path)
    ft, _ = FileType.objects.get_or_create(postfix=postfix)

    # File.path is unique -> upsert by path
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
        logger.info(f"Attached NEW File {file_path} to sample {sample} (checksum={checksum}, postfix={postfix})")
        return

    # Existing file row: update if needed (checksum may change if file regenerated)
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
        logger.info(f"Updated existing File {file_path} for sample {sample} (checksum={checksum}, postfix={postfix})")
    else:
        logger.debug(f"No changes for existing File {file_path} (already attached to sample {sample})")


# ---------------------------------------------------------------------
# Command
# ---------------------------------------------------------------------
class Command(BaseCommand):
    help = "Attach coinfection major/minor consensus files (FASTA/EMBL) to Sample via File model."

    def add_arguments(self, parser):
        parser.add_argument(
            "--dry-run",
            action="store_true",
            default=False,
            help="Do not write to DB, only log what would be attached.",
        )

    def handle(self, *args, **options):
        dry_run: bool = options["dry_run"]

        logger.info(f"Starting attach_coinfection_consensus_files (dry_run={dry_run})")
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

            major_path, minor_path = pick_major_minor_files(folder)

            if major_path is None and minor_path is None:
                logger.warning(f"Sample {sample}: no *_consensus_major.* or *_consensus_minor.* found in {folder}")
                continue

            if major_path is None:
                logger.warning(f"Sample {sample}: MAJOR consensus missing in {folder}; will attach only MINOR if present")
            if minor_path is None:
                logger.warning(f"Sample {sample}: MINOR consensus missing in {folder}; will attach only MAJOR if present")

            with transaction.atomic():
                if major_path is not None:
                    upsert_file(sample=sample, file_path=major_path, dry_run=dry_run)
                if minor_path is not None:
                    upsert_file(sample=sample, file_path=minor_path, dry_run=dry_run)

        logger.info("Done.")