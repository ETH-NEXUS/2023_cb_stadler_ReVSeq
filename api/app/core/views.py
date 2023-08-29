from django_filters import rest_framework as filters
from .models import SampleCount, Plate, Substrain, Metadata, Sample
from .serializers import (
    SampleCountSerializers,
    PlateSerializer,
    MetadataSerializer,
    SubstrainSerializer,
    AggregatedCountSerializer,
    SampleSerializer,
)
from rest_framework import viewsets
from rest_framework import filters as drf_filters
from rest_framework_csv.renderers import CSVRenderer
from rest_framework.settings import api_settings
from rest_framework.renderers import JSONRenderer
from rest_framework.decorators import action
from rest_framework.response import Response
from easydict import EasyDict as edict
from rest_framework.exceptions import ValidationError
from copy import deepcopy
from rest_framework import status


def def_value():
    return {}


class SampleCountViewSet(viewsets.ModelViewSet):
    queryset = SampleCount.objects.all()
    serializer_class = SampleCountSerializers
    filter_backends = (
        filters.DjangoFilterBackend,
        drf_filters.OrderingFilter,
    )
    filterset_fields = (
        "plate__barcode",
        "substrain__name",
        "substrain__strain__name",
        "sample__pseudoanonymized_id",
    )
    renderer_classes = (
        JSONRenderer,
        CSVRenderer,
    ) + tuple(api_settings.DEFAULT_RENDERER_CLASSES)

    # to get csv with curl:
    # curl -H "Accept: text/csv" http://localhost:8000/api/samplecounts/ > sample_counts.csv

    @action(detail=False, methods=["get"])
    def aggregate(self, request):
        pseudoanonymized_id = request.query_params.get(
            "sample__pseudoanonymized_id", None
        )
        if pseudoanonymized_id is None:
            raise ValidationError({"error": "pseudoanonymized_id filter is required"})

        try:
            sample = Sample.objects.get(pseudoanonymized_id=pseudoanonymized_id)
        except Sample.DoesNotExist:
            raise ValidationError(
                {"error": "Sample with this pseudoanonymized_id does not exist"}
            )

        queryset = self.filter_queryset(self.get_queryset())
        response_data = {
            "sample": {
                "sample_id": sample.pseudoanonymized_id,
                "plate": sample.plate.barcode,
            },
            "strains": [],
        }
        new_data = edict(
            {
                "aligned": 0,
                "length": 0,
                "rpkm": 0,
                "rpkm_proportions": 0,
                "normcounts": 0,
                "outlier": False,
                "qc_status": "",
                "coverage_threshold": 0,
                "coverage": 0,
            }
        )
        strains = {}
        for item in queryset:
            strain = item.substrain.strain.name
            if strain not in strains:
                strains[strain] = deepcopy(new_data)
            strains[strain].aligned += item.aligned
            strains[strain].length += item.length
            strains[strain].rpkm += item.rpkm
            strains[strain].rpkm_proportions += item.rpkm_proportions
            strains[strain].normcounts += item.normcounts
            strains[strain].outlier = item.outlier
            # strains[strain].qc_status = item.qc_status
            # strains[strain].qc_threshold = item.qc_threshold
            # strains[strain].coverage += item.coverage

        response_data["strains"] = [
            {"strain": key, **value} for key, value in strains.items()
        ]

        serializer = AggregatedCountSerializer(data=response_data)

        if serializer.is_valid():
            return Response(serializer.data)
        else:
            return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def finalize_response(self, request, response, *args, **kwargs):
        if isinstance(request.accepted_renderer, CSVRenderer):
            filename = "sample_counts.csv"
            response["Content-Disposition"] = f'attachment; filename="{filename}"'
        return super().finalize_response(request, response, *args, **kwargs)


class MetadataViewSet(viewsets.ModelViewSet):
    queryset = Metadata.objects.all()
    serializer_class = MetadataSerializer
    filter_backends = (filters.DjangoFilterBackend, drf_filters.OrderingFilter)
    filterset_fields = ("plate__barcode",)

    renderer_classes = (
        JSONRenderer,
        CSVRenderer,
    ) + tuple(api_settings.DEFAULT_RENDERER_CLASSES)

    # to get csv with curl:
    # curl -H "Accept: text/csv" http://localhost:8000/api/metadata/ > metadata.csv

    def finalize_response(self, request, response, *args, **kwargs):
        if isinstance(request.accepted_renderer, CSVRenderer):
            filename = "metadata.csv"
            response["Content-Disposition"] = f'attachment; filename="{filename}"'
        return super().finalize_response(request, response, *args, **kwargs)


class PlateViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = Plate.objects.all()
    serializer_class = PlateSerializer


class SubstrainViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = Substrain.objects.all()
    serializer_class = SubstrainSerializer
