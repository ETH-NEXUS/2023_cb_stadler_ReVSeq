from __future__ import annotations

import hashlib
from pathlib import Path

from django.core.management.base import BaseCommand
from django.db import transaction

from helpers.color_log import logger
from core.models import Sample, File, FileType


# ---------------------------------------------------------------------
# Hardcoded input data: files that must be RE-assigned with NEW checksum
# (delete existing for that sample+kind, then create fresh)
# ---------------------------------------------------------------------
EMBL_REASSIGNMENTS: list[dict[str, str]] = [
    {
        "sample_id": "m2-9eyPAv",
        "path": "/data/RVSeqPlate10-m2/sample_m2-9eyPAv/m2-9eyPAv_consensus_major.embl.gz",
        "kind": "major",
    },
    {
        "sample_id": "m2-F2HifQ",
        "path": "/data/RVSeqPlate5-m2/sample_m2-F2HifQ/m2-F2HifQ_consensus_major.embl.gz",
        "kind": "major",
    },
    {
        "sample_id": "m2-SrqwxQ",
        "path": "/data/RVSeqPlate11-m2/sample_m2-SrqwxQ/m2-SrqwxQ_consensus_major.embl.gz",
        "kind": "major",
    },
    {
        "sample_id": "m2-w9v3Qv",
        "path": "/data/RVSeqPlate5-m2/sample_m2-w9v3Qv/m2-w9v3Qv_consensus_minor.embl.gz",
        "kind": "minor",
    },
    {
        "sample_id": "m2-EJ7fEz",
        "path": "/data/RVSeqPlate5-m2/sample_m2-EJ7fEz/m2-EJ7fEz_consensus_minor.embl.gz",
        "kind": "minor",
    },
    {
        "sample_id": "m2-4L3pXZ",
        "path": "/data/RVSeqPlate6-m2/sample_m2-4L3pXZ/m2-4L3pXZ_consensus_major.embl.gz",
        "kind": "major",
    },
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
    lowered = filename.lower()
    if lowered.endswith(".embl.gz"):
        return ".embl.gz"
    # fallback
    return Path(filename).suffix or lowered


def is_target_embl_file(path_str: str, kind: str) -> bool:
    """
    Determines whether an existing File row is the "same kind" (major/minor)
    EMBL file we want to replace for this sample.
    """
    name = Path(path_str).name
    if not name.endswith(".embl.gz"):
        return False
    if kind == "major":
        return "_consensus_major.embl.gz" in name
    if kind == "minor":
        return "_consensus_minor.embl.gz" in name
    return False


def delete_existing_kind(*, sample: Sample, kind: str, dry_run: bool) -> int:
    """
    Delete any existing EMBL file rows for this sample that match the requested kind.
    Returns number of rows deleted.
    """
    qs = File.objects.filter(sample=sample)
    to_delete_ids: list[int] = [
        f.id for f in qs.only("id", "path") if is_target_embl_file(f.path, kind)
    ]

    if not to_delete_ids:
        return 0

    if dry_run:
        logger.info(
            f"[DRY RUN] Would delete {len(to_delete_ids)} existing {kind.upper()} EMBL File(s) for sample {sample}"
        )
        return len(to_delete_ids)

    deleted, _ = File.objects.filter(id__in=to_delete_ids).delete()
    return deleted


def create_fresh_file(*, sample: Sample, file_path: Path, dry_run: bool) -> None:
    postfix = detect_postfix(file_path.name)
    checksum = md5sum(file_path)

    if dry_run:
        logger.info(
            f"[DRY RUN] Would create NEW File row for {file_path} on sample {sample} "
            f"(checksum={checksum}, postfix={postfix})"
        )
        return

    ft, _ = FileType.objects.get_or_create(postfix=postfix)

    # If path is unique and exists globally, we must delete that row first,
    # otherwise create() will fail. This can happen if the file was previously
    # attached to the wrong sample/plate.
    existing_by_path = File.objects.filter(path=str(file_path)).first()
    if existing_by_path is not None:
        logger.warning(
            f"File path already existed in DB (path={file_path}); deleting it before re-creating."
        )
        existing_by_path.delete()

    File.objects.create(
        path=str(file_path),
        checksum=checksum,
        type=ft,
        sample=sample,
        plate=sample.plate,
    )

    logger.info(
        f"Created NEW File {file_path} for sample {sample} "
        f"(checksum={checksum}, postfix={postfix})"
    )


# ---------------------------------------------------------------------
# Command
# ---------------------------------------------------------------------
class Command(BaseCommand):
    help = (
        "Re-assign specific consensus *.embl.gz files: "
        "delete existing major/minor EMBL file(s) for the sample, then create fresh with new checksum."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "--dry-run",
            action="store_true",
            default=False,
            help="Do not write to DB, only log what would happen.",
        )

    def handle(self, *args, **options):
        dry_run: bool = options["dry_run"]

        logger.info(f"Starting reassign_embl_consensus_files (dry_run={dry_run})")
        logger.info(f"Processing {len(EMBL_REASSIGNMENTS)} EMBL reassignment entries")

        for item in EMBL_REASSIGNMENTS:
            sample_id = item["sample_id"]
            path_str = item["path"]
            kind = item["kind"]  # "major" or "minor"

            file_path = Path(path_str)

            sample = (
                Sample.objects.filter(pseudonymized_id=sample_id)
                .select_related("plate")
                .first()
            )
            if sample is None:
                logger.warning(f"Sample not found in DB: {sample_id} (path={file_path})")
                continue

            if not file_path.exists() or not file_path.is_file():
                logger.warning(f"File missing or not a file: {file_path} (sample={sample})")
                continue

            if not file_path.name.endswith(".embl.gz"):
                logger.warning(f"Not an .embl.gz file, skipping: {file_path} (sample={sample})")
                continue

            with transaction.atomic():
                deleted = delete_existing_kind(sample=sample, kind=kind, dry_run=dry_run)
                if deleted:
                    logger.info(
                        f"{sample} ({sample_id}): deleted {deleted} existing {kind.upper()} EMBL File(s) before re-create"
                    )
                else:
                    logger.info(f"{sample} ({sample_id}): no existing {kind.upper()} EMBL File found to delete")

                create_fresh_file(sample=sample, file_path=file_path, dry_run=dry_run)

        logger.info("Done.")