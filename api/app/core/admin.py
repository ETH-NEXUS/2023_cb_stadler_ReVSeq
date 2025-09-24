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
CDSPositions,
CDSCount
)
from unfold.admin import ModelAdmin

# ============== Inlines =================


class FileInline(admin.TabularInline):
    model = File
    fk_name = "sample"
    extra = 0

class SampleCountInline(admin.TabularInline):
    model = SampleCount
    fk_name = "sample"
    extra = 0

class CDSCountInline(admin.TabularInline):
    model = CDSCount
    fk_name = "sample"
    extra = 0

class MetadataInline(admin.TabularInline):
    model = Metadata
    fk_name = "sample"
    extra = 0

# Register your models here.

@admin.register(CDSCount)
class CDSCountAdmin(ModelAdmin):
    list_display = ("sample", "CDS_name", "CDS_name")
    search_fields = ("substrain", "sample", "CDS_name" )

@admin.register(CDSPositions)
class CDSPositionsAdmin(ModelAdmin):
    list_display = ("gen_bank_id", "cds_start", "cds_end")
    search_fields = ("gen_bank_id", )



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
    list_display = ("pseudonymized_id", 'valid', "sample_number", "well", "job_id", "control", "upload_to_ena")
    search_fields = ("sample_number", "pseudonymized_id", "well__location")
    list_filter_submit = True
    list_filter = ("plate__barcode", "control", "control_type")
    inlines =  [FileInline,  MetadataInline]
    filter_horizontal = ('secondary_strains',)


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
    list_display = ("name", "strain", "taxon_id", "serotype")

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
        "substrain",

    )
    list_filter_submit = True
    list_filter = ("plate__barcode",  "substrain__name")
    search_fields = ( "substrain__name", 'sample__pseudonymized_id', "sample__pseudonymized_id",)


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
