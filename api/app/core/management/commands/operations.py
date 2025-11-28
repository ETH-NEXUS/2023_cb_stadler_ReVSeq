from pathlib import Path

from django.core.management.base import BaseCommand, CommandError

from core.models import Sample
from helpers.color_log import logger


class Command(BaseCommand):
    help = (
        "Reset job_id and analysis_job_id to None for samples listed in a text file.\n"
        "The file should contain one pseudonymized_id per line."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "file_path",
            type=str,
            help="Path to the text file with pseudonymized_ids (one per line).",
        )

    def handle(self, *args, **options):
        file_path = options["file_path"]
        path = Path(file_path)

        if not path.exists():
            raise CommandError(f"File not found: {path}")

        # Read IDs, strip whitespace, ignore empty lines
        ids = [
            line.strip()
            for line in path.read_text().splitlines()
            if line.strip()
        ]

        if not ids:
            self.stdout.write(self.style.WARNING("No IDs found in file. Nothing to do."))
            return

        self.stdout.write(f"Read {len(ids)} pseudonymized IDs from {path}.")

        updated = 0
        not_found = []

        for pid in ids:
            qs = Sample.objects.filter(pseudonymized_id=pid)
            if not qs.exists():
                not_found.append(pid)
                logger.warning(f"No Sample found for pseudonymized_id={pid}")
                continue

            count = qs.update(job_id=None, analysis_job_id=None)
            updated += count
            logger.info(f"Reset job_id and analysis_job_id for {count} sample(s) with pseudonymized_id={pid}")

        self.stdout.write(self.style.SUCCESS(f"Updated {updated} Sample objects."))

        if not_found:
            self.stdout.write(
                self.style.WARNING(
                    f"{len(not_found)} IDs had no matching Sample:\n" +
                    "\n".join(not_found)
                )
            )