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


class SampleCountViewSet(viewsets.ModelViewSet):
    queryset = SampleCount.objects.all()
    serializer_class = SampleCountSerializers

    def get_queryset(self):
        queryset = SampleCount.objects.all()
        barcode = self.request.query_params.get("barcode", None)
        substrain = self.request.query_params.get("substrain", None)
        if barcode is not None:
            queryset = queryset.filter(plate__barcode=barcode)
        if substrain is not None:
            queryset = queryset.filter(substrain__name=substrain)
        return queryset


class MetadataViewSet(viewsets.ModelViewSet):
    queryset = Metadata.objects.all()
    serializer_class = MetadataSerializer

    def get_queryset(self):
        queryset = Metadata.objects.all()
        barcode = self.request.query_params.get("barcode", None)
        if barcode is not None:
            queryset = queryset.filter(plate__barcode=barcode)
        return queryset


class PlateListView(generics.ListAPIView):
    queryset = Plate.objects.all()
    serializer_class = PlateSerializer


class SubstrainListView(generics.ListAPIView):
    queryset = Substrain.objects.all()
    serializer_class = SubstrainSerializer
