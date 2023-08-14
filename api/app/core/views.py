from django.shortcuts import render
from django_filters import rest_framework as filters
from .models import SampleCount, Plate, Substrain, Metadata
from .serializers import (
    SampleCountSerializers,
    PlateSerializer,
    MetadataSerializer,
    SubstrainSerializer,
)
from rest_framework import viewsets
from rest_framework import generics
from rest_framework import filters as drf_filters


class SampleCountViewSet(viewsets.ModelViewSet):
    queryset = SampleCount.objects.all()
    serializer_class = SampleCountSerializers
    filter_backends = (filters.DjangoFilterBackend, drf_filters.OrderingFilter)
    filterset_fields = ("plate__barcode", "substrain__name", "substrain__strain__name")


class MetadataViewSet(viewsets.ModelViewSet):
    queryset = Metadata.objects.all()
    serializer_class = MetadataSerializer
    filter_backends = (filters.DjangoFilterBackend, drf_filters.OrderingFilter)
    filterset_fields = ("plate__barcode",)

    # def get_queryset(self):
    #     queryset = Metadata.objects.all()
    #     barcode = self.request.query_params.get("barcode", None)
    #     if barcode is not None:
    #         queryset = queryset.filter(plate__barcode=barcode)
    #     return queryset


class PlateListView(generics.ListAPIView):
    queryset = Plate.objects.all()
    serializer_class = PlateSerializer


class SubstrainListView(generics.ListAPIView):
    queryset = Substrain.objects.all()
    serializer_class = SubstrainSerializer
