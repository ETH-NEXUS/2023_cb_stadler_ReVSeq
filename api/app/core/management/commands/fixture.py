from django.core.management.base import BaseCommand
import traceback
import csv
from core.models import Panel, Strain, Substrain, FileType

FILE_PANEL_STRAIN = "initial_data/panel_strain.csv"
FILE_STRAIN_SUBSTRAIN = "initial_data/strain_substrain.csv"
FILE_FILE_TYPE = "initial_data/file_type.csv"


def read_csv_file(input_file):
    with open(input_file, "r", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "what", type=str, help="Which fixture to upload"
        )  # panel_strain, strain_substrain, file_type

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

    def handle(self, *args, **options):
        try:
            if options.get("what") == "panel_strain":
                self.panel_strain()
            elif options.get("what") == "strain_substrain":
                self.strain_substrain()
            elif options.get("what") == "file_type":
                self.file_type()
        except Exception as ex:
            print(ex)
            traceback.print_exc()
