import hashlib
from pathlib import Path

from django.core.management.base import BaseCommand
from django.db import transaction

from core.models import Sample, File, FileType
from helpers.color_log import logger


class Command(BaseCommand):
    help = "Associate filtered EMBL files with Samples (and fix wrong attachments if needed)"

    # ----------------------------------------------
    # 1) static mapping (sample_id → embl_file_path)
    # ----------------------------------------------
    FILE_MAP = [
        ("m2-33tHNX", "/data/RVSeqPlate13-m2/sample_m2-33tHNX/m2-33tHNX_filtered.embl.gz"),
        ("m2-3MgVni", "/data/RVSeqPlate3-m2/sample_m2-3MgVni/m2-3MgVni_filtered.embl.gz"),
        ("m2-3URpmc", "/data/RVSeqPlate3-m2/sample_m2-3URpmc/m2-3URpmc_filtered.embl.gz"),
        ("m2-4Rt9sw", "/data/RVSeqPlate3-m2/sample_m2-4Rt9sw/m2-4Rt9sw_filtered.embl.gz"),
        ("m2-BdtP8F", "/data/RVSeqPlate3-m2/sample_m2-BdtP8F/m2-BdtP8F_filtered.embl.gz"),
        ("m2-Qb4JVo", "/data/RVSeqPlate5-m2/sample_m2-Qb4JVo/m2-Qb4JVo_filtered.embl.gz"),
        ("m2-xxAxC9", "/data/RVSeqPlate5-m2/sample_m2-xxAxC9/m2-xxAxC9_filtered.embl.gz"),
        ("m2-voFSf6", "/data/RVSeqPlate6-m2/sample_m2-voFSf6/m2-voFSf6_filtered.embl.gz"),
        ("m2-os283P", "/data/RVSeqPlate7-m2/sample_m2-os283P/m2-os283P_filtered.embl.gz"),
        ("m2-XMdsys", "/data/RVSeqPlate7-m2/sample_m2-XMdsys/m2-XMdsys_filtered.embl.gz"),
        ("m2-Y4wzR3", "/data/RVSeqPlate7-m2/sample_m2-Y4wzR3/m2-Y4wzR3_filtered.embl.gz"),
        ("m2-4Q8Stt", "/data/RVSeqPlate8-m2/sample_m2-4Q8Stt/m2-4Q8Stt_filtered.embl.gz"),
        ("m2-BbET5x", "/data/RVSeqPlate8-m2/sample_m2-BbET5x/m2-BbET5x_filtered.embl.gz"),
        ("m2-LYceeY", "/data/RVSeqPlate8-m2/sample_m2-LYceeY/m2-LYceeY_filtered.embl.gz"),
        ("m2-NESNdH", "/data/RVSeqPlate8-m2/sample_m2-NESNdH/m2-NESNdH_filtered.embl.gz"),
    ]

    # ----------------------------------------------------
    # MD5 checksum helper
    # ----------------------------------------------------
    def compute_md5(self, file_path: str) -> str:
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()

    # ----------------------------------------------------
    # Main execution
    # ----------------------------------------------------
    @transaction.atomic
    def handle(self, *args, **options):
        logger.info("Starting EMBL file → Sample association / fix operation")

        # ensure FileType exists
        postfix = ".embl.gz"
        filetype, _ = FileType.objects.get_or_create(postfix=postfix)

        processed = 0
        missing_samples = []
        missing_files = []
        created_files = []
        reassigned_files = []

        for sample_id, file_path in self.FILE_MAP:
            # check sample exists
            sample = Sample.objects.filter(pseudonymized_id=sample_id).first()
            if not sample:
                missing_samples.append(sample_id)
                logger.error(f"Sample {sample_id} NOT FOUND")
                continue

            # check file exists on disk
            if not Path(file_path).exists():
                missing_files.append(file_path)
                logger.error(f"File does not exist: {file_path}")
                continue

            # compute checksum
            checksum = self.compute_md5(file_path)

            # get or create File by path
            file_obj, created = File.objects.get_or_create(
                path=file_path,
                defaults={
                    "checksum": checksum,
                    "type": filetype,
                    "sample": sample,
                    "plate": sample.plate,
                },
            )

            if created:
                logger.info(f"Created new File for sample {sample_id}: {file_path}")
                created_files.append(file_path)
            else:
                # already existed – possibly linked to wrong sample
                changed = False

                if file_obj.sample != sample:
                    logger.warning(
                        "Reassigning File %s from sample %s to sample %s",
                        file_path,
                        file_obj.sample.pseudonymized_id if file_obj.sample else "None",
                        sample.pseudonymized_id,
                    )
                    file_obj.sample = sample
                    changed = True

                if file_obj.plate != sample.plate:
                    logger.info(
                        "Updating plate for File %s from %s to %s",
                        file_path,
                        file_obj.plate.barcode if file_obj.plate else "None",
                        sample.plate.barcode if sample.plate else "None",
                    )
                    file_obj.plate = sample.plate
                    changed = True

                if file_obj.type != filetype:
                    logger.info(
                        "Updating FileType for %s from %s to %s",
                        file_path,
                        file_obj.type.postfix if file_obj.type else "None",
                        filetype.postfix,
                    )
                    file_obj.type = filetype
                    changed = True

                # keep checksum in sync with actual file on disk
                if file_obj.checksum != checksum:
                    logger.info(
                        "Updating checksum for %s (old=%s, new=%s)",
                        file_path,
                        file_obj.checksum,
                        checksum,
                    )
                    file_obj.checksum = checksum
                    changed = True

                if changed:
                    file_obj.save()
                    reassigned_files.append(file_path)
                else:
                    logger.info(f"File entry already correct: {file_path}")

            processed += 1

        logger.info("\n---- SUMMARY ----")
        logger.info(f"Processed entries: {processed}")
        logger.info(f"New File objects created: {len(created_files)}")
        logger.info(f"Reassigned/updated File objects: {len(reassigned_files)}")

        if missing_samples:
            logger.warning(f"Missing samples: {missing_samples}")
        if missing_files:
            logger.warning(f"Missing files: {missing_files}")

        logger.info("Done.")