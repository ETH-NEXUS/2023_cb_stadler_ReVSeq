from .models import (
    Plate,
    Well,
    Sample,
    Strain,
    Panel,
    Substrain,
    SampleCount,
    Metadata,
    File,
    FileType,
CDSPositions
)

from rest_framework import serializers

class CDSPositionsSerializer(serializers.ModelSerializer):
    class Meta:
        model = CDSPositions
        fields = "__all__"


class FileSerializer(serializers.ModelSerializer):
    class Meta:
        model = File
        fields = (
            "id",
            "path",
            "checksum",
        )


class PlateSerializer(serializers.ModelSerializer):
    files = FileSerializer(many=True, read_only=True, source="file_set")

    class Meta:
        model = Plate
        fields = ("id", "barcode", "files")


class WellSerializer(serializers.ModelSerializer):
    plate = PlateSerializer(read_only=True)

    class Meta:
        model = Well
        fields = ("location", "plate")


class SampleSerializer(serializers.ModelSerializer):
    files = FileSerializer(many=True, read_only=True, source="file_set")
    well = WellSerializer(read_only=True)
    plate = PlateSerializer(read_only=True)

    class Meta:
        model = Sample
        fields = ("well", "pseudonymized_id", "plate", "files", "control", "control_type")


class StrainSerializer(serializers.ModelSerializer):
    class Meta:
        model = Strain
        fields = "__all__"


class PanelSerializer(serializers.ModelSerializer):
    strain = StrainSerializer(read_only=True)

    class Meta:
        model = Panel
        fields = ("name", "strain")


class SubstrainSerializer(serializers.ModelSerializer):
    strain = StrainSerializer(read_only=True)

    class Meta:
        model = Substrain
        fields = ("name", "strain")


class SampleCountSerializers(serializers.ModelSerializer):
    plate = PlateSerializer(read_only=True)
    substrain = SubstrainSerializer(read_only=True)
    sample = SampleSerializer(read_only=True)
    strain = serializers.SerializerMethodField("get_strain")

    def get_strain(self, obj):
        return obj.substrain.strain.name

    class Meta:
        model = SampleCount
        fields = (
            "plate",
            "sample",
            "strain",
            "substrain",
            "aligned",
            "length",
            "rpkm",
            "rpkm_proportions",
            "normcounts",
            "outlier",
            "coverage_threshold",
            "coverage",
            "coverage_status",
            "readnum_status",
            "readnum_threshold",
            "percentile_threshold",
            "tax_id",
            "scientific_name",
            "DP",
            "consensus_number_n",
            "consensus_fraction_n",
            "consensus",
            "consensus_cds",

        )


class MetadataSerializer(serializers.ModelSerializer):
    plate = PlateSerializer(read_only=True)
    sample = SampleSerializer(read_only=True)
    well = WellSerializer(read_only=True)

    class Meta:
        model = Metadata
        fields = (
            "plate",
            "sample",
            "well",
            "prescriber",
            "order_date",
            "ent_date",
            "treatment_type",
            "data",
        )


"""
 ---------------- For aggregation -----------------
"""


class StrainCountsSerializer(serializers.Serializer):
    aligned = serializers.IntegerField()
    length = serializers.IntegerField()
    rpkm = serializers.FloatField()
    rpkm_proportions = serializers.FloatField()
    normcounts = serializers.FloatField()
    outlier = serializers.BooleanField()
    strain = serializers.CharField()
    coverage_threshold = serializers.FloatField()
    coverage = serializers.FloatField()
    coverage_status = serializers.CharField()
    readnum_status = serializers.CharField()
    readnum_threshold = serializers.FloatField()
    percentile_threshold = serializers.CharField()
    tax_id = serializers.IntegerField()
    scientific_name = serializers.CharField()
    DP20 = serializers.CharField()




class SimpleSampleSerializer(serializers.Serializer):
    sample_id = serializers.CharField()
    plate = serializers.CharField()

    class Meta:
        model = Sample
        fields = "pseudonymized_id"


class AggregatedCountSerializer(serializers.Serializer):
    sample = SimpleSampleSerializer(many=False)
    strains = StrainCountsSerializer(many=True)

    class Meta:
        fields = ("sample", "strains")
