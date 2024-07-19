from django.core.handlers.wsgi import WSGIRequest
from django.views import View
from django_filters import rest_framework as filters
from drf_spectacular.utils import extend_schema
from rest_framework.views import APIView

from .models import SampleCount, Plate, Substrain, Metadata, Sample, File, CDSPositions
from .serializers import (
    SampleCountSerializers,
    PlateSerializer,
    MetadataSerializer,
    SubstrainSerializer,
    AggregatedCountSerializer,
    SampleSerializer,
    FileSerializer,
CDSPositionsSerializer
)
from rest_framework import viewsets
from rest_framework import filters as drf_filters
from rest_framework_csv.renderers import CSVRenderer
from rest_framework.settings import api_settings
from rest_framework.decorators import action
from rest_framework.response import Response
from easydict import EasyDict as edict
from rest_framework.exceptions import ValidationError
from copy import deepcopy
from rest_framework import status
from django.http import FileResponse, Http404
import os
from rest_framework.decorators import api_view

from rest_framework.renderers import JSONRenderer
from rest_framework.permissions import IsAuthenticated
from rest_framework_simplejwt.authentication import JWTAuthentication
from rest_framework.authentication import SessionAuthentication
from rest_framework.decorators import api_view, authentication_classes, permission_classes
import logging
import json
from django.http import JsonResponse
from django.core.management import call_command
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt
from rest_framework.request import Request


logger = logging.getLogger(__name__)




def def_value():
    return {}


class SampleViewSet(viewsets.ModelViewSet):
    """
    API Endpoint: Samples

    This endpoint provides access to the samples within a plate. Each sample
    is associated with a specific well on a plate and contains an array of file paths.
    These files can be downloaded using the designated download endpoint.

    Features:
    - List View: Retrieve a list of all available samples.
    - Detail View: Access details of a specific sample by providing
      its unique database ID.

    Filters:
    - By Plate Barcode: Filter samples based on the barcode of the plate they are associated with.
    - By pseudonymized ID: Locate a specific sample using its pseudonymized ID.

    Allowed HTTP Methods: GET, HEAD, OPTIONS
    """
    authentication_classes = [JWTAuthentication, SessionAuthentication]
    permission_classes = [IsAuthenticated]
    queryset = Sample.objects.all()
    serializer_class = SampleSerializer
    filter_backends = (
        filters.DjangoFilterBackend,
        drf_filters.OrderingFilter,
    )
    filterset_fields = ("plate__barcode", "pseudonymized_id", "control")

    http_method_names = ["get", "head", "options"]

    @action(detail=False, methods=["post"])
    def ena_upload_study(self, request, *args, **kwargs):
        """
        API Endpoint: ENA Upload Study

        This endpoint initiates the upload of the study data to the ENA (European Nucleotide Archive).

        Features:
        - Study Upload: Initiates the upload of the study data to the ENA.

        Allowed HTTP Methods: POST
        """
        call_command('ena_upload', 'study')
        return JsonResponse({"detail": "Study upload initiated."}, status=200)

    @action(detail=False, methods=["post"])
    def ena_upload_ser_and_analysis(self, request, *args, **kwargs):

        """
        API Endpoint: ENA Upload SER and Analysis

        This endpoint initiates the upload of the SER and Analysis data to the ENA (European Nucleotide Archive).

        Features:
        - SER and Analysis Upload: Initiates the upload of the SER and Analysis data to the ENA.

        Allowed HTTP Methods: POST
        """
        call_command('ena_upload', 'ser_and_analysis')
        return JsonResponse({"detail": "SER and Analysis upload initiated."}, status=200)


class SampleCountViewSet(viewsets.ModelViewSet):
    """
    API Endpoint: Sample Counts

    This endpoint provides a view into the counts associated with each sample.
    Sample counts  indicate various attributes such as alignment, length,
    rpkm (Reads Per Kilobase Million), and more.

    Features:
    - List View: Access a list of all available sample counts.
    - Aggregate View: Provides aggregated metrics for a given sample.
     This action requires the pseudonymized_id  of the sample to be provided in the query parameter  to fetch relevant metrics ( 'sample__pseudonymized_id') .

    Filters:
    - By Plate Barcode: Filter counts associated with a specific plate's barcode.
    - By Substrain Name: View counts specific to a particular substrain.
    - By Strain Name: Filter based on the overarching strain.
    - By pseudonymized_id: Filter counts using the sample's pseudonymized_id.

    Supported Formats:
    - JSON

    CSV File Retrieval:
    - **Swagger Caveat**: Within the Swagger interface, selecting the CSV format will yield the error "Could not satisfy the request Accept header."
        However, the CSV data can be obtained using a designated header via `curl`.
            Example:
            ```bash
            curl -H "Accept: text/csv" http://https://revseq.nexus.ethz.ch/api/samplecounts/ > sample_counts.csv
            ```
    - Alternatively, a csv file can be downloaded directly from our default DRF API view.

    Allowed HTTP Methods: GET, HEAD, OPTIONS

    """
    authentication_classes = [JWTAuthentication, SessionAuthentication]
    permission_classes = [IsAuthenticated]

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
        "sample__pseudonymized_id",
    )
    renderer_classes = (
        JSONRenderer,
        CSVRenderer,
    ) + tuple(api_settings.DEFAULT_RENDERER_CLASSES)

    http_method_names = ["get", "head", "options"]

    def get_queryset(self):
        return SampleCount.objects.prefetch_related('plate', 'substrain', 'sample')

    """
    Call the aggregate function it like this, with the obligatory parameter sample__pseudonymized_id:
    
     /api/samplecounts/aggregate/?sample__pseudonymized_id=896a06
    """

    @action(
        detail=False,
        methods=["get"],
        url_path="aggregate/(?P<sample__pseudonymized_id>[^/.]+)",
    )
    def aggregate(self, request, sample__pseudonymized_id=None):

        if not sample__pseudonymized_id:
            raise ValidationError({"error": "pseudonymized_id is required"})

        try:
            sample = Sample.objects.get(pseudonymized_id=sample__pseudonymized_id)

        except Sample.DoesNotExist:
            raise ValidationError(
                {"error": "Sample with this pseudonymized_id does not exist"}
            )

        queryset = SampleCount.objects.filter(sample=sample)   # SampleCount.objects.filter(sample=sample)

        response_data = {
            "sample": {
                "sample_id": sample.pseudonymized_id,
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
                "coverage_threshold": 0,
                "coverage": 0,
                "readnum_threshold": 0,




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
            strains[strain].coverage_threshold = item.coverage_threshold
            strains[strain].coverage += item.coverage
            strains[strain].readnum_threshold +=item.readnum_threshold



        response_data["strains"] = [
            {"strain": key, **value} for key, value in strains.items()
        ]
        #print("RESPONSE DATA", response_data)

        serializer = AggregatedCountSerializer(data=response_data)
        print("SERIALIZER", serializer)

        if serializer.is_valid():
            return Response(serializer.data)
        else:
            return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def finalize_response(self, request, response, *args, **kwargs):
        if hasattr(request, "accepted_renderer") and isinstance(
            request.accepted_renderer, CSVRenderer
        ):
            filename = "sample_counts.csv"
            response["Content-Disposition"] = f'attachment; filename="{filename}"'
        return super().finalize_response(request, response, *args, **kwargs)


class MetadataViewSet(viewsets.ModelViewSet):
    """
    API Endpoint: Metadata

    This endpoint provides access to metadata related to the samples or plates.

    Features:
    - List View: Access a list of all available metadata entries.
    - Filters:
      - By Plate Barcode: Filter metadata associated with a specific plate's barcode.

    Supported Formats:
    - JSON

    CSV File Retrieval:
    - **Swagger Caveat**: Within the Swagger interface, selecting the CSV format will yield the error
      "Could not satisfy the request Accept header."
      However, the CSV data can be obtained using a designated header via `curl`.
          Example:
          ```bash
          curl -H "Accept: text/csv" http://https://revseq.nexus.ethz.ch//api/metadata/ > metadata.csv
          ```
      - Alternatively, a CSV file can be downloaded directly  from our default DRF API view.

    Allowed HTTP Methods: GET, HEAD, OPTIONS

    """
    authentication_classes = [JWTAuthentication, SessionAuthentication]
    permission_classes = [IsAuthenticated]
    queryset = Metadata.objects.all()
    serializer_class = MetadataSerializer
    filter_backends = (filters.DjangoFilterBackend, drf_filters.OrderingFilter)
    filterset_fields = ("plate__barcode",)

    renderer_classes = (
        JSONRenderer,
        CSVRenderer,
    ) + tuple(api_settings.DEFAULT_RENDERER_CLASSES)
    http_method_names = ["get", "head", "options"]

    def finalize_response(self, request, response, *args, **kwargs):
        if isinstance(request.accepted_renderer, CSVRenderer):
            filename = "metadata.csv"
            response["Content-Disposition"] = f'attachment; filename="{filename}"'
        return super().finalize_response(request, response, *args, **kwargs)


from django.core.exceptions import ObjectDoesNotExist


class PlateViewSet(viewsets.ModelViewSet):
    """
    API Endpoint: Plates

    This endpoint provides access to plates' data.

    Features:
    - List View: Retrieve a list of all available plates.
    - Filtered View: Obtain details of a specific plate by providing its barcode.

    Filters:
    - By Plate Barcode: Obtain details of a specific plate using its barcode.
      Example:
      ```
      /api/plates/?barcode=PLATE_BARCODE_HERE
      ```
      If a plate with the provided barcode exists, details of that plate will be returned.
      If not, an empty list is returned.

    Supported Formats:
    - JSON

    Allowed HTTP Methods: GET, HEAD, OPTIONS

    Note: The barcode is not a primary key in this viewset, hence the custom filtering.
    """

    authentication_classes = [JWTAuthentication, SessionAuthentication]
    permission_classes = [IsAuthenticated]
    queryset = Plate.objects.all()
    serializer_class = PlateSerializer
    filter_backends = (
        filters.DjangoFilterBackend,
        drf_filters.OrderingFilter,
    )
    filterset_fields = ("barcode",)
    http_method_names = ["get", "head", "options"]


class SubstrainViewSet(viewsets.ReadOnlyModelViewSet):
    """
    API Endpoint: Substrains

    This endpoint provides access to the data of different substrains.

    Features:
    - List View: Retrieve a comprehensive list of all available substrains in the system.
    - Detail View: Access detailed information of a specific substrain by using its unique ID.

    Supported Formats:
    - JSON

    Allowed HTTP Methods: GET, HEAD, OPTIONS
    """

    authentication_classes = [JWTAuthentication, SessionAuthentication]
    permission_classes = [IsAuthenticated]

    queryset = Substrain.objects.all()
    serializer_class = SubstrainSerializer
    http_method_names = ["get", "head", "options"]


class CDSPositionsViewSet(viewsets.ModelViewSet):
    """
    API Endpoint: CDS Positions

    """
    authentication_classes = [JWTAuthentication, SessionAuthentication]
    permission_classes = [IsAuthenticated]
    queryset = CDSPositions.objects.all()
    serializer_class = CDSPositionsSerializer
    http_method_names = ["get", "head", "options"]

    filter_backends = (
        filters.DjangoFilterBackend,
        drf_filters.OrderingFilter,
    )

    filterset_fields = (
        "sample__pseudonymized_id",
    )




class FileViewSet(viewsets.ModelViewSet):
    """
    API Endpoint: Files

    Overview:
    The endpoint facilitates access to file paths tied to a specific sample or plate.
    The files can be downloaded using the designated download endpoint.



    Available Filters:
    - **Plate Barcode**: Narrow down file paths associated with a particular plate's barcode.
    - **Sample pseudonymized ID**: Filter files specific to a given sample using its unique pseudonymized ID.

    Supported Data Formats:
    - JSON

    Permitted HTTP Methods: GET, HEAD, OPTIONS"""

    authentication_classes = [JWTAuthentication, SessionAuthentication]
    permission_classes = [IsAuthenticated]

    queryset = File.objects.all()
    serializer_class = FileSerializer
    http_method_names = ["get", "head", "options"]

    filter_backends = (
        filters.DjangoFilterBackend,
        drf_filters.OrderingFilter,
    )

    filterset_fields = (
        "plate__barcode",
        "sample__pseudonymized_id",
    )

@api_view(["GET"])
@authentication_classes([JWTAuthentication, SessionAuthentication])
@permission_classes([IsAuthenticated])
def download_file(request, filepath):
    """
    API Endpoint: Download File

    This endpoint offers a direct download mechanism for a file based on the provided filepath.

    Features:
    - Direct File Download: By supplying a valid filepath as a parameter, users can initiate the download of the specified file.

    Parameters:
    - filepath: Absolute path of the desired file (a list of available files is provided with each plate or sample item, as well as within the files API endpoint.)

    Response:
    - If the file exists and is accessible: Initiates a file download.
    - If the file does not exist or is inaccessible: Returns an appropriate error message.

    Allowed HTTP Methods: GET
    """

    file_path = filepath
    if not file_path.startswith("/"):
        file_path = "/" + file_path

    if os.path.exists(file_path):
        try:
            response = FileResponse(
                open(file_path, "rb"), content_type="application/octet-stream"
            )
            response[
                "Content-Disposition"
            ] = f"attachment; filename={os.path.basename(file_path)}"
            return response
        except PermissionError:
            return Response(
                {"detail": "Permission denied for the requested file."},
                status=status.HTTP_403_FORBIDDEN,
            )
        except FileNotFoundError:
            return Response(
                {"detail": "File not found."},
                status=status.HTTP_404_NOT_FOUND,
            )
    else:
        raise Http404("File doesn't exist or is inaccessible.")

@method_decorator(csrf_exempt, name='dispatch')
@extend_schema(exclude=True)
class ImportResultsView(APIView):
    """
    post:
    Import results data of the experiment.

    path -- The file system path to the directory containing the data to import.
    """
    def post(self, request: Request, *args, **kwargs):
        logger.debug("Starting import")
        path = request.data.get("path")
        if not path:
            return Response({"detail": "Please provide path."}, status=status.HTTP_400_BAD_REQUEST)

        try:
            call_command('import', "--import_dir", path)
            logger.debug("Import finished")
            return Response({"detail": "Import finished."}, status=status.HTTP_200_OK)
        except Exception as e:
            logger.error("Import failed: %s", e, exc_info=True)
            return Response({"detail": f"Import failed: {e}"}, status=status.HTTP_400_BAD_REQUEST)