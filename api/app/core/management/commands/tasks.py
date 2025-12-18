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


# root@api:/app# ./manage.py tasks
# API:  INFO ----- Starting attach_coinfection_consensus_files (dry_run=False)
# API:  INFO ----- Processing 23 sample folder entries
# API:  INFO ----- Attached NEW File /data/RVSeqPlate10-m2/sample_m2-Dq37e8/m2-Dq37e8_consensus_major.fa.gz to sample m2-Dq37e8 (checksum=4829af56704ac484f7cee107c12baaf4, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate10-m2/sample_m2-Dq37e8/m2-Dq37e8_consensus_minor.fa.gz to sample m2-Dq37e8 (checksum=60c8abc520365095e4da8d7045439538, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate10-m2/sample_m2-wfEdma/m2-wfEdma_consensus_major.fa.gz to sample m2-wfEdma (checksum=4d5de43fc3862c8587b38d82f2abe4c4, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate10-m2/sample_m2-wfEdma/m2-wfEdma_consensus_minor.fa.gz to sample m2-wfEdma (checksum=f0ff6b719a1334d958188270f3a51466, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate10-m2/sample_m2-9eyPAv/m2-9eyPAv_consensus_major.embl.gz to sample m2-9eyPAv (checksum=8045bc6f80c3e5b18f2208e5304492a1, postfix=.embl.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate10-m2/sample_m2-9eyPAv/m2-9eyPAv_consensus_minor.fa.gz to sample m2-9eyPAv (checksum=da05873d07f321b030a832b7ca360bcc, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate11-m2/sample_m2-GRjNzo/m2-GRjNzo_consensus_major.fa.gz to sample m2-GRjNzo (checksum=eb40763e9f5d5545798924c0f5717d4a, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate11-m2/sample_m2-GRjNzo/m2-GRjNzo_consensus_minor.fa.gz to sample m2-GRjNzo (checksum=773aefd207ddc018107433ca468dcdca, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate11-m2/sample_m2-Kg3Fdm/m2-Kg3Fdm_consensus_major.fa.gz to sample m2-Kg3Fdm (checksum=5ce986432ce16e314b1bc5789d787071, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate11-m2/sample_m2-Kg3Fdm/m2-Kg3Fdm_consensus_minor.fa.gz to sample m2-Kg3Fdm (checksum=144bd13608c259fc5c4579d19fdabb0e, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate11-m2/sample_m2-SrqwxQ/m2-SrqwxQ_consensus_major.embl.gz to sample m2-SrqwxQ (checksum=c35dc50a34dbb000f4e840e3d667daa1, postfix=.embl.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate11-m2/sample_m2-SrqwxQ/m2-SrqwxQ_consensus_minor.fa.gz to sample m2-SrqwxQ (checksum=92bd5a182395dd52d487bbfe89b99103, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate13-m2/sample_m2-mzMRtV/m2-mzMRtV_consensus_major.fa.gz to sample m2-mzMRtV (checksum=919b26276d00ab848c2063a9887843e5, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate13-m2/sample_m2-mzMRtV/m2-mzMRtV_consensus_minor.fa.gz to sample m2-mzMRtV (checksum=3e9fd1394f07e320b347c7b957d5ab26, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate3-m2/sample_m2-zwqknf/m2-zwqknf_consensus_major.fa.gz to sample m2-zwqknf (checksum=b1b7d2fa1adeed8a87df30853a208dbf, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate3-m2/sample_m2-zwqknf/m2-zwqknf_consensus_minor.fa.gz to sample m2-zwqknf (checksum=271e5039bd2daa51297a46ca721299fa, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate3-m2/sample_m2-3oZJ5x/m2-3oZJ5x_consensus_major.fa.gz to sample m2-3oZJ5x (checksum=bdabe52023241de36bba2579df9ca07f, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate3-m2/sample_m2-3oZJ5x/m2-3oZJ5x_consensus_minor.fa.gz to sample m2-3oZJ5x (checksum=0a0b760c3c2e68808ed7fc33ac8ca722, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate3-m2/sample_m2-43Vmvj/m2-43Vmvj_consensus_major.fa.gz to sample m2-43Vmvj (checksum=1f53c8d4d30a10145b974e417ef48b18, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate3-m2/sample_m2-43Vmvj/m2-43Vmvj_consensus_minor.fa.gz to sample m2-43Vmvj (checksum=65707f17778bf13d018fe3d7a8536668, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate4-m2/sample_m2-NmHGDa/m2-NmHGDa_consensus_major.fa.gz to sample m2-NmHGDa (checksum=c6576d5a03e7e49afee17a2780c27dba, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate4-m2/sample_m2-NmHGDa/m2-NmHGDa_consensus_minor.fa.gz to sample m2-NmHGDa (checksum=0488eba9ebef4d04c4324ecc0ef1984f, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate4-m2/sample_m2-u89dBs/m2-u89dBs_consensus_major.fa.gz to sample m2-u89dBs (checksum=7767f76da41f0adfc74c8d0439818cd1, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate4-m2/sample_m2-u89dBs/m2-u89dBs_consensus_minor.fa.gz to sample m2-u89dBs (checksum=e354a513c32cf496fdf9045420828cb8, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate4-m2/sample_m2-gnDemF/m2-gnDemF_consensus_major.fa.gz to sample m2-gnDemF (checksum=6a0404b24d1e7dcde86e0d22ee6829f2, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate4-m2/sample_m2-gnDemF/m2-gnDemF_consensus_minor.fa.gz to sample m2-gnDemF (checksum=6777fce16b201be3633f3e688b32f09d, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate4-m2/sample_m2-HXHQAT/m2-HXHQAT_consensus_major.fa.gz to sample m2-HXHQAT (checksum=15818ec6ad5fc9f9c4b28d78cd7971b0, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate4-m2/sample_m2-HXHQAT/m2-HXHQAT_consensus_minor.fa.gz to sample m2-HXHQAT (checksum=281a5a93982af8eae78292bc60840b21, postfix=.fa.gz)
# API:  WARNING ----- Sample not found in DB: m2-RVSeqPlate5-KOpos (folder=/data/RVSeqPlate5-m2/sample_m2-RVSeqPlate5-KOpos)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate5-m2/sample_m2-w9v3Qv/m2-w9v3Qv_consensus_major.fa.gz to sample m2-w9v3Qv (checksum=1d2f451d0b77b0de208ecafb8880b94d, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate5-m2/sample_m2-w9v3Qv/m2-w9v3Qv_consensus_minor.embl.gz to sample m2-w9v3Qv (checksum=8b104cfcce905f10dc8715234526d30a, postfix=.embl.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate5-m2/sample_m2-F2HifQ/m2-F2HifQ_consensus_major.embl.gz to sample m2-F2HifQ (checksum=ee1ac1f598fd1a0c8c93af4c241946f5, postfix=.embl.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate5-m2/sample_m2-F2HifQ/m2-F2HifQ_consensus_minor.fa.gz to sample m2-F2HifQ (checksum=48bfc28c8964d6a0de6603dbd96fa632, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate5-m2/sample_m2-EJ7fEz/m2-EJ7fEz_consensus_major.fa.gz to sample m2-EJ7fEz (checksum=734e9a3356bd9d47ad5834ba73a57d9b, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate5-m2/sample_m2-EJ7fEz/m2-EJ7fEz_consensus_minor.embl.gz to sample m2-EJ7fEz (checksum=0473d5cb08816027e8f781716583e436, postfix=.embl.gz)
# API:  WARNING ----- Sample not found in DB: m2-RVSeqPlate6-KOpos (folder=/data/RVSeqPlate6-m2/sample_m2-RVSeqPlate6-KOpos)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate6-m2/sample_m2-4L3pXZ/m2-4L3pXZ_consensus_major.embl.gz to sample m2-4L3pXZ (checksum=754aa7b663db06bb7ff796c01c12d971, postfix=.embl.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate6-m2/sample_m2-4L3pXZ/m2-4L3pXZ_consensus_minor.fa.gz to sample m2-4L3pXZ (checksum=f3006125010fe13f2a59d9a53e5321c7, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate7-m2/sample_m2-NDEjHn/m2-NDEjHn_consensus_major.fa.gz to sample m2-NDEjHn (checksum=53f11672fb9bf3986d810a0a8e349e57, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate7-m2/sample_m2-NDEjHn/m2-NDEjHn_consensus_minor.fa.gz to sample m2-NDEjHn (checksum=63054e9b83d8a5a80ffc839163de1d37, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate8-m2/sample_m2-Y4rZ7q/m2-Y4rZ7q_consensus_major.fa.gz to sample m2-Y4rZ7q (checksum=9c22352c302a68301cc606e0fbee7e49, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate8-m2/sample_m2-Y4rZ7q/m2-Y4rZ7q_consensus_minor.fa.gz to sample m2-Y4rZ7q (checksum=6448b77ed1a7da54b7dc587de3eab7ce, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate8-m2/sample_m2-XHdiUs/m2-XHdiUs_consensus_major.fa.gz to sample m2-XHdiUs (checksum=6b0fbaee34ede98c1670416aef4af5c6, postfix=.fa.gz)
# API:  INFO ----- Attached NEW File /data/RVSeqPlate8-m2/sample_m2-XHdiUs/m2-XHdiUs_consensus_minor.fa.gz to sample m2-XHdiUs (checksum=db3c35174faca484a661bab57f0e8ae8, postfix=.fa.gz)