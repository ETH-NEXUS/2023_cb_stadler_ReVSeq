from django.contrib import admin
from .models import (
    Plate,
    Well,
    Sample,
    SampleFile,
    Strain,
    Panel,
    Substrain,
    SampleCount,
    Metadata,
    FileType,
    File,
)

# Register your models here.


@admin.register(Plate)
class PlateAdmin(admin.ModelAdmin):
    list_display = ("barcode",)
    search_fields = ("barcode",)


@admin.register(Well)
class WellAdmin(admin.ModelAdmin):
    list_display = ("location", "plate")
    search_fields = ("location",)


@admin.register(Sample)
class SampleAdmin(admin.ModelAdmin):
    list_display = ("pseudoanonymized_id", "sample_number", "well")
    search_fields = ("sample_number",)


@admin.register(SampleFile)
class SampleFileAdmin(admin.ModelAdmin):
    list_display = ("sample", "file_type", "file")
    search_fields = ("sample__sample_id", "file_type", "file")


@admin.register(Strain)
class StrainAdmin(admin.ModelAdmin):
    list_display = ("name",)
    search_fields = ("name",)


@admin.register(Panel)
class PanelAdmin(admin.ModelAdmin):
    list_display = ("name", "strain")
    search_fields = ("name", "strain__name")


@admin.register(Substrain)
class SubstrainAdmin(admin.ModelAdmin):
    list_display = ("name", "strain")

    search_fields = ("name", "strain__name")


@admin.register(SampleCount)
class SampleCountAdmin(admin.ModelAdmin):
    list_display = (
        "aligned",
        "length",
        "rpkm",
        "rpkm_proportions",
        "normcounts",
        "outlier",
    )


@admin.register(Metadata)
class MetadataAdmin(admin.ModelAdmin):
    list_display = (
        "well",
        "plate",
        "sample",
        "prescriber",
    )


@admin.register(FileType)
class FileTypeAdmin(admin.ModelAdmin):
    list_display = ("postfix",)


@admin.register(File)
class FileAdmin(admin.ModelAdmin):
    list_display = ("path", "checksum", "type", "sample", "plate")
    search_fields = (
        "path",
        "checksum",
        "type__postfix",
        "sample__sample_number",
        "plate__barcode",
    )
