from django.core.management.base import BaseCommand
from django.db import transaction
from core.models import File


OLD_PREFIX = "/data/revseq/results/gather_results/"
NEW_PREFIX = "/data/"
SUFFIX = "consensus_upload.gz"


class Command(BaseCommand):
    help = "Fix file paths by replacing old prefix with /data/. Use --dry-run to only print changes."

    def add_arguments(self, parser):
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Only print what would be changed, do not save anything.",
        )

    def handle(self, *args, **options):
        dry_run = options["dry_run"]

        qs = File.objects.filter(
            path__startswith=OLD_PREFIX,
            path__endswith=SUFFIX,
        )

        if not qs.exists():
            self.stdout.write(self.style.WARNING("No matching file paths found."))
            return

        self.stdout.write(f"Found {qs.count()} files to process...")

        # Don't wrap in transaction for dry-run
        context = transaction.atomic() if not dry_run else None

        if context:
            context.__enter__()

        try:
            for f in qs:
                old_path = f.path
                new_path = old_path.replace(OLD_PREFIX, NEW_PREFIX, 1)

                if dry_run:
                    self.stdout.write(
                        f"[DRY RUN] {old_path}\n          → {new_path}\n"
                    )
                else:
                    f.path = new_path
                    f.save(update_fields=["path"])
                    self.stdout.write(f"Updated:\n  {old_path}\n  → {new_path}")

        finally:
            # Roll back automatically for dry run
            if dry_run and context:
                context.__exit__(None, None, None)

        if dry_run:
            self.stdout.write(self.style.WARNING("Dry run complete — no changes saved."))
        else:
            self.stdout.write(self.style.SUCCESS("File path updates saved."))