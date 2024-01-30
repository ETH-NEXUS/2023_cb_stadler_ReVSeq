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
from unfold.admin import ModelAdmin


# Register your models here.




@admin.register(Plate)
class PlateAdmin(ModelAdmin): # class PlateAdmin(admin.ModelAdmin):
    list_display = ("barcode",)
    search_fields = ("barcode",)


@admin.register(Well)
class WellAdmin(ModelAdmin):
    list_display = ("location", "plate")
    search_fields = ("location",)


@admin.register(Sample)
class SampleAdmin(ModelAdmin):
    list_display = ("pseudonymized_id", "sample_number", "well", "job_id")
    search_fields = ("sample_number", "pseudonymized_id", "well__location")
    list_filter_submit = True
    list_filter = ("plate__barcode",)


@admin.register(SampleFile)
class SampleFileAdmin(ModelAdmin):
    list_display = ("sample", "file_type", "file")
    search_fields = ("sample__sample_id", "file_type", "file")


@admin.register(Strain)
class StrainAdmin(ModelAdmin):
    list_display = ("name",)
    search_fields = ("name",)


@admin.register(Panel)
class PanelAdmin(ModelAdmin):
    list_display = ("name", "strain")
    search_fields = ("name", "strain__name")


@admin.register(Substrain)
class SubstrainAdmin(ModelAdmin):
    list_display = ("name", "strain")

    search_fields = ("name", "strain__name")
    list_filter_submit = True
    list_filter= ("strain",)


@admin.register(SampleCount)
class SampleCountAdmin(ModelAdmin):
    list_display = (
        "aligned",
        "length",
        "rpkm",
        "rpkm_proportions",
        "normcounts",
        "outlier",
    )
    list_filter_submit = True
    list_filter = ("plate__barcode", "sample__pseudonymized_id", "substrain__name")


@admin.register(Metadata)
class MetadataAdmin(ModelAdmin):
    list_display = (
        "well",
        "plate",
        "sample",
        "prescriber",
    )
    list_filter_submit = True
    list_filter = ("plate__barcode", )


@admin.register(FileType)
class FileTypeAdmin(ModelAdmin):
    list_display = ("postfix",)


@admin.register(File)
class FileAdmin(ModelAdmin):
    list_display = ("path", "checksum", "type", "sample", "plate")
    search_fields = (
        "path",
        "checksum",
        "type__postfix",
        "sample__sample_number",
        "plate__barcode",
    )
