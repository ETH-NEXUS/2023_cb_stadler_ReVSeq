from django.core.management.base import BaseCommand
import traceback
import csv
from core.models import Panel, Strain, Substrain

FILE_PANEL_STRAIN_LOOKUP = "data/panel_strain_lookup.csv"
FILE_STRAIN_SUBSTRAIN_LOOKUP = "data/strain_substrain_lookup.csv"


def read_csv_file(input_file):
    with open(input_file, "r", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "what", type=str, help="Which fixture to upload"
        )  # panel_strain, strain_substrain

    def strain_substrain(self):
        data = read_csv_file(FILE_STRAIN_SUBSTRAIN_LOOKUP)
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
        data = read_csv_file(FILE_PANEL_STRAIN_LOOKUP)
        for item in data:
            if item["panel_name"]:
                panel, _ = Panel.objects.get_or_create(name=item["panel_name"])
                strain, _ = Strain.objects.get_or_create(name=item["strain_name"])
                panel.strain = strain
                panel.save()
                strain.save()

    def handle(self, *args, **options):
        try:
            if options.get("what") == "panel_strain":
                self.panel_strain()
            elif options.get("what") == "strain_substrain":
                self.strain_substrain()
        except Exception as ex:
            print(ex)
            traceback.print_exc()
