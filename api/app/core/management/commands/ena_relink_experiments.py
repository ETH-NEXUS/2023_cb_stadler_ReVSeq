import os
import time
import xml.etree.ElementTree as ET
from typing import List
import requests
from django.core.management.base import BaseCommand, CommandError


# docker compose exec api python manage.py ena_relink_experiments \
#   --webin-user "user" \
#   --webin-pass "password" \
#   --old-study ERP167225 \
#   --new-study ERP192961 \
#   --new-study-secondary PRJEB112404 \
#   --dry-run


class Command(BaseCommand):
    help = "One-off command: relink existing ENA experiments to the new raw read study."

    def add_arguments(self, parser):
        parser.add_argument("--webin-user", required=True)
        parser.add_argument("--webin-pass", required=True)
        parser.add_argument("--old-study", default="ERP167225")
        parser.add_argument("--new-study", default="ERP192961")
        parser.add_argument("--new-study-secondary", required=True)
        parser.add_argument(

            "--ignore-sample",
            action="append",
            default=[],
            help="Sample accession to ignore. Can be used multiple times.",

        )
        parser.add_argument(
            "--dev",
            action="store_true",
            help="Use ENA dev endpoints instead of production.",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Create XML but do not submit it.",
        )

    def handle(self, *args, **options):
        webin_user = options["webin_user"]
        webin_pass = options["webin_pass"]
        old_study = options["old_study"]
        new_study = options["new_study"]
        new_study_secondary = options["new_study_secondary"]
        ignored_samples = set(options["ignore_sample"])
        use_dev = options["dev"]
        dry_run = options["dry_run"]
        if use_dev:
            ena_report_url = "https://wwwdev.ebi.ac.uk/ena/submit/report"
            ena_dropbox_url = "https://wwwdev.ebi.ac.uk/ena/submit/drop-box"
        else:
            ena_report_url = "https://www.ebi.ac.uk/ena/submit/report"
            ena_dropbox_url = "https://www.ebi.ac.uk/ena/submit/drop-box"
        self.stdout.write(f"Old study: {old_study}")
        self.stdout.write(f"New raw read study: {new_study}")
        self.stdout.write(f"New raw read secondary accession: {new_study_secondary}")
        self.stdout.write(f"Using {'DEV' if use_dev else 'PRODUCTION'} ENA endpoints")

        experiment_ids = self.fetch_experiment_ids(
            ena_report_url=ena_report_url,
            webin_user=webin_user,
            webin_pass=webin_pass,
            old_study=old_study,
        )

        self.stdout.write(f"Found {len(experiment_ids)} experiments")
        if not experiment_ids:
            raise CommandError("No experiments found. Stopping.")
        experiment_xml_strings = []
        for experiment_id in experiment_ids:
            self.stdout.write(f"Processing {experiment_id}...")
            current_xml = self.fetch_experiment_xml(
                ena_dropbox_url=ena_dropbox_url,
                webin_user=webin_user,
                webin_pass=webin_pass,
                experiment_id=experiment_id,
            )
            sample_accession = self.extract_sample_accession(current_xml)
            if sample_accession in ignored_samples:
                self.stdout.write(f"Skipping ignored sample {sample_accession}")
                continue
            modified_xml = self.modify_experiment_study_ref(
                current_xml=current_xml,
                new_study=new_study,
                new_study_secondary=new_study_secondary,
            )

            experiment_xml_strings.append(modified_xml)
            # Be gentle with ENA.
            time.sleep(0.2)
        if not experiment_xml_strings:
            raise CommandError("No experiments left after filtering ignored samples.")
        experiments_xml = self.build_experiment_set_xml(experiment_xml_strings)
        submission_xml = self.build_modify_submission_xml()
        self.write_file("experiments.xml", experiments_xml)
        self.write_file("modify.xml", submission_xml)
        self.stdout.write("Created experiments.xml")
        self.stdout.write("Created modify.xml")
        if dry_run:
            self.stdout.write(self.style.WARNING("Dry run enabled. Not submitting."))
            return
        receipt_text = self.submit_modified_experiments(
            ena_dropbox_url=ena_dropbox_url,
            webin_user=webin_user,
            webin_pass=webin_pass,
            submission_xml_path="modify.xml",
            experiments_xml_path="experiments.xml",

        )
        self.write_file("modify_experiments_receipt.xml", receipt_text)
        self.stdout.write("Created modify_experiments_receipt.xml")
        if 'success="true"' in receipt_text:
            self.stdout.write(self.style.SUCCESS("ENA modify submission succeeded."))
        else:
            self.stdout.write(self.style.ERROR("ENA modify submission may have failed. Check receipt."))

    def fetch_experiment_ids(

            self,

            *,

            ena_report_url: str,

            webin_user: str,

            webin_pass: str,

            old_study: str,

    ) -> List[str]:

        """
        Fetch all public experiment IDs from the old study.
        ENA may return only a limited number of results per request,
        so we ask for a larger page size here.
        """
        url = f"{ena_report_url}/experiments"
        params = {
            "study_accession": old_study,
            "status": "PUBLIC",
            "format": "json",
            "maxresults": "5000",

        }
        response = requests.get(
            url,
            params=params,
            auth=(webin_user, webin_pass),
            timeout=60,
        )
        response.raise_for_status()
        data = response.json()
        experiment_ids = []
        for item in data:
            report = item.get("report", {})
            experiment_id = report.get("id")
            if experiment_id:
                experiment_ids.append(experiment_id)
        return experiment_ids

    def fetch_experiment_xml(
            self,
            *,
            ena_dropbox_url: str,
            webin_user: str,
            webin_pass: str,
            experiment_id: str,

    ) -> str:

        url = f"{ena_dropbox_url}/experiments/{experiment_id}"
        response = requests.get(
            url,
            auth=(webin_user, webin_pass),
            timeout=60,
        )
        response.raise_for_status()
        return response.text

    def extract_sample_accession(self, current_xml: str) -> str | None:

        root = ET.fromstring(current_xml)
        sample_descriptor = root.find(".//SAMPLE_DESCRIPTOR")
        if sample_descriptor is None:
            return None
        return sample_descriptor.attrib.get("accession")

    def modify_experiment_study_ref(
            self,
            *,
            current_xml: str,
            new_study: str,
            new_study_secondary: str,
    ) -> str:
        root = ET.fromstring(current_xml)
        # The fetched XML:
        # <EXPERIMENT_SET><EXPERIMENT>...</EXPERIMENT></EXPERIMENT_SET>
        # or directly <EXPERIMENT>...</EXPERIMENT>.
        if root.tag == "EXPERIMENT_SET":
            experiment = root.find("EXPERIMENT")
        elif root.tag == "EXPERIMENT":
            experiment = root
        else:
            raise CommandError(f"Unexpected root element: {root.tag}")
        if experiment is None:
            raise CommandError("Could not find EXPERIMENT element.")
        old_study_ref = experiment.find("STUDY_REF")
        if old_study_ref is None:
            raise CommandError("Could not find STUDY_REF element.")
        old_study_ref.clear()
        old_study_ref.attrib["accession"] = new_study
        identifiers = ET.SubElement(old_study_ref, "IDENTIFIERS")
        primary_id = ET.SubElement(identifiers, "PRIMARY_ID")
        primary_id.text = new_study
        secondary_id = ET.SubElement(identifiers, "SECONDARY_ID")
        secondary_id.text = new_study_secondary
        return ET.tostring(experiment, encoding="unicode")

    def build_experiment_set_xml(self, experiment_xml_strings: List[str]) -> str:
        lines = ["<EXPERIMENT_SET>"]
        for experiment_xml in experiment_xml_strings:
            lines.append(experiment_xml)
        lines.append("</EXPERIMENT_SET>")
        return "\n".join(lines)

    def build_modify_submission_xml(self) -> str:

        return """<SUBMISSION>
  <ACTIONS>
    <ACTION>
      <MODIFY/>
    </ACTION>
  </ACTIONS>
</SUBMISSION>

"""

    def submit_modified_experiments(
            self,
            *,
            ena_dropbox_url: str,
            webin_user: str,
            webin_pass: str,
            submission_xml_path: str,
            experiments_xml_path: str,
    ) -> str:
        url = f"{ena_dropbox_url}/submit/"
        with open(submission_xml_path, "rb") as submission_file:
            with open(experiments_xml_path, "rb") as experiments_file:
                files = {
                    "SUBMISSION": submission_file,
                    "EXPERIMENT": experiments_file,
                }
                response = requests.post(
                    url,
                    files=files,
                    auth=(webin_user, webin_pass),
                    timeout=120,
                )
        response.raise_for_status()
        return response.text

    def write_file(self, path: str, content: str) -> None:
        with open(path, "w", encoding="utf-8") as file:
            file.write(content)
