
import sys
from django.core.management.base import BaseCommand
from helpers.color_log import logger
from core.importer import Importer

class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument(
            "--import_dir", "-d", type=str, help="Directory to import from"
        )

    def handle(self, *args, **options):
        if len(sys.argv) < 2:
            logger.info("Usage: python import.py -d <import_dir>")
            sys.exit(1)
        import_dir = options.get("import_dir")
        importer = Importer()
        importer.handle_import(import_dir)


