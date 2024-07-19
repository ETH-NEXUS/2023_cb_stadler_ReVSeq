from django.core.management.base import BaseCommand
import traceback
import csv
from core.models import Panel, Strain, Substrain, FileType, CDSPositions
from helpers.color_log import logger

FILE_PANEL_STRAIN = "initial_data/panel_strain.csv"
FILE_STRAIN_SUBSTRAIN = "initial_data/strain_substrain.csv"
FILE_FILE_TYPE = "initial_data/file_type.csv"
FILE_TAXON_ID = "initial_data/substrain_taxon_lookup.csv"
FILE_CDS = "initial_data/cds.bed"


def read_csv_file(input_file):
    with open(input_file, "r", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "what", type=str, help="Which fixture to upload"
        )  # panel_strain, strain_substrain, file_type, taxon_id, cds

    def strain_substrain(self):
        data = read_csv_file(FILE_STRAIN_SUBSTRAIN)
        for item in data:
            if item["strain_name"]:
                strain, _ = Strain.objects.get_or_create(name=item["strain_name"])
                substrain, _ = Substrain.objects.get_or_create(
                    name=item["substrain_name"]
                )
                substrain.strain = strain
                substrain.save()
                strain.save()

    def substrain_taxon_id(self):
        data = read_csv_file(FILE_TAXON_ID)
        for item in data:
            if item["substrain_name"]:
                substrain, _ = Substrain.objects.get_or_create(
                    name=item["substrain_name"]
                )
                substrain.taxon_id = item["tax_id"]
                substrain.scientific_name = item["scientific_name"]
                substrain.save()
                print(f"Updated {substrain}")


    def panel_strain(self):
        data = read_csv_file(FILE_PANEL_STRAIN)
        for item in data:
            if item["panel_name"]:
                panel, _ = Panel.objects.get_or_create(name=item["panel_name"])
                strain, _ = Strain.objects.get_or_create(name=item["strain_name"])
                panel.strain = strain
                panel.save()
                strain.save()

    def file_type(self):
        data = read_csv_file(FILE_FILE_TYPE)
        for item in data:
            postfix = item["postfix"]
            file_type, _ = FileType.objects.get_or_create(postfix=postfix)
            print(f"Created {file_type}")

    def cds(self):
        try:
            with open(FILE_CDS, "r") as file:
                data = file.readlines()

            for item in data:
                item_list = item.strip().split("\t")
                if len(item_list) < 3:
                    logger.error(f"Invalid line format: {item}")
                    continue

                try:
                    gen_bank_id = item_list[0].strip()
                    cds_start = int(item_list[1].strip())
                    cds_end = int(item_list[2].strip())

                    defaults = {
                        "cds_start": cds_start,
                        "cds_end": cds_end,
                    }

                    cds_position, created = CDSPositions.objects.update_or_create(
                        gen_bank_id=gen_bank_id, defaults=defaults
                    )
                    logger.info(f"{'Created' if created else 'Updated'} CDS position: {cds_position}")

                except ValueError as ve:
                    logger.error(f"Value error for line {item}: {ve}")
                except Exception as e:
                    logger.error(f"Error processing line {item}: {e}")

        except Exception as e:
            logger.error(f"Error reading file {FILE_CDS}: {e}")

    def handle(self, *args, **options):
        try:
            if options.get("what") == "panel_strain":
                self.panel_strain()
            elif options.get("what") == "strain_substrain":
                self.strain_substrain()
            elif options.get("what") == "file_type":
                self.file_type()
            elif options.get("what") == "taxon_id":
                self.substrain_taxon_id()
            elif options.get("what") == "cds":
                self.cds()
        except Exception as ex:
            print(ex)
            traceback.print_exc()
