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
)

from rest_framework import serializers


class PlateSerializer(serializers.ModelSerializer):
    class Meta:
        model = Plate
        fields = "__all__"


class WellSerializer(serializers.ModelSerializer):
    plate = PlateSerializer(read_only=True)

    class Meta:
        model = Well
        fields = ("location", "plate")


class FileTypeSerializer(serializers.ModelSerializer):
    class Meta:
        model = FileType
        fields = "__all__"


class FileSerializer(serializers.ModelSerializer):
    type = FileTypeSerializer(read_only=True)

    class Meta:
        model = File
        fields = ("path", "checksum", "type")


class SampleSerializer(serializers.ModelSerializer):
    files = FileSerializer(many=True, read_only=True, source="file_set")
    well = WellSerializer(read_only=True)
    plate = PlateSerializer(read_only=True)

    class Meta:
        model = Sample
        fields = ("well", "pseudoanonymized_id", "plate", "files")


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
            "qc_status",
            "coverage_threshold",
            "coverage",
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
    qc_status = serializers.CharField()
    coverage_threshold = serializers.FloatField()
    coverage = serializers.FloatField()


class SimpleSampleSerializer(serializers.Serializer):
    sample_id = serializers.CharField()
    plate = serializers.CharField()

    class Meta:
        model = Sample
        fields = "pseudoanonymized_id"


class AggregatedCountSerializer(serializers.Serializer):
    sample = SimpleSampleSerializer(many=False)
    strains = StrainCountsSerializer(many=True)

    class Meta:
        fields = ("sample", "strains")
