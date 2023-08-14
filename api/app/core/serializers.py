from .models import (
    Plate,
    Well,
    Sample,
    Strain,
    Panel,
    Substrain,
    SampleCount,
    Metadata,
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


class SampleSerializer(serializers.ModelSerializer):
    well = WellSerializer(read_only=True)

    class Meta:
        model = Sample
        fields = ("well", "pseudoanonymized_id")


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
