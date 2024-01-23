from rest_framework.test import APIClient, APITestCase
from core.models import Sample, SampleCount
from django.urls import reverse
from django.core.management import call_command
import os

dir_name = os.path.dirname(os.path.abspath(__file__))

metadata_file_path = os.path.join(dir_name, "test_data", "010623RVS_metadata.csv")
counts_file_path = os.path.join(dir_name, "test_data", "0b9ad2_count_table.csv")
URL_SAMPLE_COUNTS_LIST = reverse("samplecounts-list")
AGGREGATE_URL = reverse("samplecounts-aggregate")


class SampleCountViewSetTest(APITestCase):
    def setUp(self):
        call_command("fixture", "panel_strain")
        call_command("fixture", "strain_substrain")
        call_command(
            "import",  # we use the very first example data for testing
            metadata_file=metadata_file_path,
            counts_file=counts_file_path,
        )
        self.client = APIClient()

    def test_get_sample_count(self):
        response = self.client.get(URL_SAMPLE_COUNTS_LIST)
        self.assertEqual(response.status_code, 200)

    def test_get_sample_count_with_filter(
        self,
    ):
        response = self.client.get(
            URL_SAMPLE_COUNTS_LIST,
            {"sample__pseudonymized_id": "0b9ad2"},
        )
        # print(response.data)
        data = response.data
        self.assertEqual(len(data), 41)
        self.assertEqual(data[0]["plate"]["barcode"], "010623RVS")
        self.assertEqual(response.status_code, 200)

    def test_aggregate_action_missing_param(self):
        response = self.client.get(AGGREGATE_URL)
        self.assertEqual(response.status_code, 400)
        self.assertIn("error", response.data)

    def test_aggregate_action(self):
        response = self.client.get(
            AGGREGATE_URL, {"sample__pseudonymized_id": "0b9ad2"}
        )
        print("#################", response.data)
