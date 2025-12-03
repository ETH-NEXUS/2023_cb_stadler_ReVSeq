"""
Management command to upload data to ENA using ENAUploader.

USAGE EXAMPLES
==============

1) Upload STUDY only
--------------------
# uses ENAUploader.upload_study()
python manage.py ena_upload --type study

2) Upload SER + ANALYSIS for specific samples (live run)
--------------------------------------------------------
python manage.py ena_upload --type ser_and_analysis --samples RyXauM 375EUk ABC123

3) Upload SER + ANALYSIS for samples from a text file
-----------------------------------------------------
samples.txt contains:
RyXauM
375EUk
ABC123

python manage.py ena_upload --type ser_and_analysis --samples-file /path/to/samples.txt

4) Upload SER ONLY (no analysis) for given samples
--------------------------------------------------
python manage.py ena_upload --type ser_and_analysis --no-analysis --samples XXXX XXXXX

5) Upload SER + ANALYSIS for INFLUENZA samples only (EMBL files only)
---------------------------------------------------------------------
python manage.py ena_upload --type ser_and_analysis --influenza-only --samples XXXX XXXX
# example for test-run:
# python manage.py ena_upload --type ser_and_analysis --influenza-only --test-run --samples m2-NESNdH

6) Upload SER (no analysis) for INFLUENZA samples only
------------------------------------------------------
python manage.py ena_upload --type ser_and_analysis --no-analysis --influenza-only --samples-file /path/to/influenza_samples.txt

7) TEST RUN - use it if you upload to the dev instance of ENA
-------------------------------------------------------------
python manage.py ena_upload --type ser_and_analysis --test-run --samples RyXauM 375EUk

8) RESEND ANALYSIS JOBS for existing SER jobs
---------------------------------------------
python manage.py ena_upload --task resend_analysis_jobs --samples RyXauM 375EUk

# Influenza-only variant:
python manage.py ena_upload --task resend_analysis_jobs --influenza-only --samples-file /path/to/samples.txt

NOTE:
- You must always specify samples (via --samples and/or --samples-file)
  for SER upload and resend_analysis_jobs.
- There is intentionally NO default "all samples with job_id is null" as it was before
"""

from pathlib import Path
from typing import List, Set

from django.core.management.base import BaseCommand, CommandError

from helpers.color_log import logger
from core.ena_uploader import ENAUploader  # adjust import path to where ENAUploader lives


class Command(BaseCommand):
    help = "Upload study/SER/analysis data to ENA or resend analysis jobs."

    def add_arguments(self, parser):
        # What to do: type vs task
        parser.add_argument(
            "-t",
            "--type",
            type=str,
            choices=["study", "ser_and_analysis"],
            help=(
                "Type of operation: "
                "'study' to upload only the study, "
                "'ser_and_analysis' to upload SER (and optionally analysis)."
            ),
        )

        parser.add_argument(
            "-r",
            "--task",
            type=str,
            choices=["resend_analysis_jobs"],
            help=(
                "Optional task to perform instead of --type.\n"
                "Currently supported: 'resend_analysis_jobs' to recreate analysis jobs "
                "for already uploaded SER jobs."
            ),
        )

        # Samples can be given directly on CLI...
        parser.add_argument(
            "-s",
            "--samples",
            nargs="*",
            help=(
                "List of pseudonymized sample IDs to process. "
                "Example: --samples RyXauM 375EUk ABC123"
            ),
        )

        # ...or via a text file (one ID per line).
        parser.add_argument(
            "--samples-file",
            type=str,
            help=(
                "Path to a text file containing one pseudonymized sample ID per line. "
                "Blank lines and lines starting with '#' are ignored."
            ),
        )

        # Flags modifying SER/analysis behavior
        parser.add_argument(
            "-a",
            "--no-analysis",
            action="store_true",
            dest="no_analysis",
            help="If set, analysis jobs will NOT be created, only SER jobs.",
        )

        parser.add_argument(
            "--influenza-only",
            action="store_true",
            dest="influenza_only",
            help=(
                "If set, analysis will only upload EMBL files (no consensus/chr). "
                "Use for influenza samples where only EMBL is required."
            ),
        )

        # test vs live run
        parser.add_argument(
            "--test-run",
            action="store_true",
            dest="test_run",
            default=False,
            help=(
                "Run in TEST mode: store job IDs in sample.test_job_id / "
                "sample.test_analysis_job_id instead of the real fields."
            ),
        )

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #

    def _load_samples_from_file(self, path_str: str) -> List[str]:
        """
        Read pseudonymized IDs line-by-line from a text file.
        Ignores empty lines and lines starting with '#'.
        """
        path = Path(path_str)
        if not path.exists():
            raise CommandError(f"Samples file does not exist: {path}")

        ids: List[str] = []
        with path.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                ids.append(line)
        return ids

    def _combine_samples(
        self,
        samples_arg: List[str] | None,
        samples_file: str | None,
    ) -> List[str]:
        """
        Combine samples from --samples and --samples-file,
        preserve order but avoid duplicates.
        """
        result: List[str] = []
        seen: Set[str] = set()

        # from CLI
        if samples_arg:
            for s in samples_arg:
                if s not in seen:
                    seen.add(s)
                    result.append(s)

        # from file
        if samples_file:
            file_ids = self._load_samples_from_file(samples_file)
            for s in file_ids:
                if s not in seen:
                    seen.add(s)
                    result.append(s)

        return result

    # ------------------------------------------------------------------ #
    # Main handler
    # ------------------------------------------------------------------ #

    def handle(self, *args, **options):
        op_type = options.get("type")
        task = options.get("task")
        samples_arg = options.get("samples")
        samples_file = options.get("samples_file")
        no_analysis = options.get("no_analysis", False)
        influenza_only = options.get("influenza_only", False)
        test_run = options.get("test_run", False)

        # Combine samples from arguments and file
        pseudonymized_ids = self._combine_samples(samples_arg, samples_file)

        if task and op_type:
            raise CommandError(
                "You cannot specify both --task and --type at the same time. "
                "Choose one."
            )

        uploader = ENAUploader(test_run=test_run)

        # --------------------- TASK: resend_analysis_jobs --------------------- #
        if task == "resend_analysis_jobs":
            if not pseudonymized_ids:
                raise CommandError(
                    "resend_analysis_jobs requires at least one sample ID "
                    "(use --samples and/or --samples-file)."
                )
            logger.info(f"Resending analysis jobs for {len(pseudonymized_ids)} samples (test_run={test_run}, influenza_only={influenza_only}).")
            uploader.resend_analysis_jobs(
                pseudonymized_ids=pseudonymized_ids,
                influenza_only=influenza_only,
            )
            return

        # ----------------------- TYPE: study / ser --------------------------- #
        if op_type == "study":
            logger.info("Uploading study to ENA (test_run=%s).", test_run)
            uploader.upload_study()
            return

        if op_type == "ser_and_analysis":
            if not pseudonymized_ids:
                raise CommandError(
                    "SER upload requires at least one sample ID "
                    "(use --samples and/or --samples-file). "
                    "There is no default filter anymore."
                )
            logger.info(f"Uploading SER for {len(pseudonymized_ids)} samples (test_run={test_run}, no_analysis={no_analysis}, influenza_only={influenza_only}).")
            uploader.upload_ser(
                pseudonymized_ids=pseudonymized_ids,
                no_analysis=no_analysis,
                influenza_only=influenza_only,
            )
            return


        raise CommandError(
            "Please specify either --type (study, ser_and_analysis) or "
            "--task (resend_analysis_jobs). See command help for examples."
        )