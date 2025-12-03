import hashlib
from pathlib import Path

from django.core.management.base import BaseCommand
from django.db import transaction

from core.models import Sample, File, FileType
from helpers.color_log import logger


SAMPLE_IDS = [
    "m2-CLwScw", "m2-4UF9K6", "m2-jtjCjd", "m2-32WNFL", "m2-375EUk",
    "m2-xwkVbK", "m2-itJLqU", "m2-g33df2", "m2-opfXjQ", "m2-Sm9GYA",
    "m2-eAr6Lo", "m2-xHGR9u", "m2-4rvG7f", "m2-SQZYHn", "m2-oLXZ9a",
    "m2-9HhvMh", "m2-ss4Fh8", "m2-VA9Yi5", "m2-qtZxkz", "m2-39kc9U",
    "m2-q2oXnd", "m2-x6NRYc", "m2-kBgpZ4", "m2-V5iihT", "m2-m6243p",
    "m2-gkcF5C", "m2-3H4uLZ", "m2-juVwzj", "m2-h2ZJd5", "m2-SpXNmN",
    "m2-K2udH5", "m2-4t7Ncr", "m2-hzNQRG", "m2-uTpHe5", "m2-BnbHLr",
    "m2-LosGHn", "m2-LCBdav", "m2-sr4bvM", "m2-iDqyV9", "m2-mmYEfJ",
    "m2-VT3bcF", "m2-sAyJww", "m2-b99ooj", "m2-2hKeW2", "m2-j3ZXpL",
    "m2-uSTYZ2", "m2-f8rfBe", "m2-FXFt6J", "m2-6z2wf2", "m2-5cJWyf",
    "m2-sFFUfJ", "m2-wVBQbT", "m2-kwbwf8", "m2-RyXauM", "m2-ydSpNQ",
    "m2-SKJgwA", "m2-DwiwwU", "m2-xYnaUy", "m2-xrRPFj", "m2-yhE88h",
    "m2-gaQkyT", "m2-TPgies"
]


class Command(BaseCommand):
    help = "Scan sample directories for *.embl.gz files and attach them as File objects."

    def compute_md5(self, file_path: Path) -> str:
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()

    @transaction.atomic
    def handle(self, *args, **options):
        logger.info("Starting EMBL file association command")

        filetype, _ = FileType.objects.get_or_create(postfix=".embl.gz")

        created_files = 0
        missing_samples = []
        no_directory_found = []
        no_embl_found = []

        for sample_id in SAMPLE_IDS:
            sample = Sample.objects.filter(pseudonymized_id=sample_id).first()
            if not sample:
                logger.error(f"❌ Sample not found: {sample_id}")
                missing_samples.append(sample_id)
                continue

            # find an existing file to derive the directory
            existing_file = File.objects.filter(sample=sample).first()
            if not existing_file:
                logger.error(f"❌ No File objects exist yet for sample {sample_id}")
                no_directory_found.append(sample_id)
                continue

            directory = Path(existing_file.path).parent

            if not directory.exists():
                logger.error(f"❌ Directory does not exist: {directory}")
                no_directory_found.append(sample_id)
                continue

            # find *.embl.gz files
            embl_files = list(directory.glob("*.embl.gz"))
            if not embl_files:
                logger.warning(f"⚠️ No .embl.gz files found for sample {sample_id} in {directory}")
                no_embl_found.append(sample_id)
                continue

            for embl_path in embl_files:
                checksum = self.compute_md5(embl_path)

                file_obj, created = File.objects.get_or_create(
                    path=str(embl_path),
                    defaults={
                        "checksum": checksum,
                        "type": filetype,
                        "sample": sample,
                        "plate": sample.plate,
                    },
                )

                if created:
                    created_files += 1
                    logger.info(f"🆕 Added EMBL file for {sample_id}: {embl_path}")
                else:
                    # ensure correct sample association
                    if file_obj.sample != sample:
                        logger.warning(f"🔄 Re-attaching file {embl_path} to sample {sample_id}")
                        file_obj.sample = sample
                        file_obj.plate = sample.plate
                        file_obj.save()

        # --- SUMMARY ---
        logger.info("\n=== SUMMARY ===")
        logger.info(f"Created new EMBL File objects: {created_files}")
        logger.info(f"Missing samples: {missing_samples}")
        logger.info(f"No directory found: {no_directory_found}")
        logger.info(f"No EMBL files found: {no_embl_found}")

        logger.info("Done.")