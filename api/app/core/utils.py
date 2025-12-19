import os
from typing import Iterable

import requests
from django.utils import timezone

from helpers.color_log import logger
from os import environ
from core.models import SampleCount, Sample, Metadata, File
import datetime as dt
# ----------------------------------------------------------------------
# HTTP helper
# ----------------------------------------------------------------------

def handle_http_request(url, payload=None, method='get', message=None, headers=None):
    try:
        response = None
        if method == 'get':
            response = requests.get(url, headers=headers)
        elif method == 'post':
            response = requests.post(url, headers=headers, json=payload)
        if response.status_code in [200, 201]:
            if message:
                logger.info(message)
            return response.json()
        else:
            logger.error(f'Error calling the endpoint {url}: {response.text}')
    except Exception as e:
        logger.error(f'Exception: {e}')
    return None


# ----------------------------------------------------------------------
# Default configuration values (read once at import time)
# ----------------------------------------------------------------------

class Defaults:
    STUDY_ENDPOINT = environ.get(
        "STUDY_ENDPOINT", "http://ena:5000/ena/api/jobs/study/"
    )
    SER_ENDPOINT = environ.get(
        "SER_ENDPOINT", "http://ena:5000/ena/api/jobs/ser/"
    )
    ANALYSIS_ENDPOINT = environ.get(
        "ANALYSIS_ENDPOINT", "http://ena:5000/ena/api/analysisjobs/"
    )
    ANALYSIS_FILES_ENDPOINT = environ.get(
        "ANALYSIS_FILES_ENDPOINT", "http://ena:5000/ena/api/analysisfiles/"
    )
    ENQUEUE_ENDPOINT = environ.get(
        "ENQUEUE_ENDPOINT", "http://ena:5000/ena/api/analysisjobs/<job_id>/enqueue/"
    )
    RELEASE_JOB_ENDPOINT = environ.get(
        "RELEASE_JOB_ENDPOINT", "http://ena:5000/ena/api/jobs/<job_id>/release/"
    )
    RELEASE_ANALYSIS_JOB_ENDPOINT = environ.get(
        "RELEASE_ANALYSIS_JOB_ENDPOINT", "http://ena:5000/ena/api/analysisjobs/<job_id>/release/"
    )
    JOBS_ENDPOINT = environ.get(
        "JOBS_ENDPOINT", "http://ena:5000/ena/api/jobs/"
    )
    ANALYSIS_JOBS_ENDPOINT = environ.get(
        "ANALYSIS_JOBS_ENDPOINT", "http://ena:5000/ena/api/analysisjobs/"
    )

    CONSENSUS_FILE_SUFFIX = environ.get(
        "CONSENSUS_FILE_SUFFIX", "consensus_upload.fasta.gz"
    )
    CHROMOSOME_MAJOR_FILE_SUFFIX = environ.get(
        "CHROMOSOME_MAJOR_FILE_SUFFIX", "chr_file_major.txt.gz"
    )
    CHROMOSOME_MINOR_FILE_SUFFIX = environ.get(
        "CHROMOSOME_MINOR_FILE_SUFFIX", "chr_file_minor.txt.gz"
    )
    CHROMOSOME_FILE_NAME = environ.get(
        "CHROMOSOME_FILE_NAME", "chr_file.txt.gz"
    )
    EMBL_FILE_SUFFIX = environ.get("EMBL_FILE_SUFFIX", "embl.gz")
    EMBL_FILE_TYPE = environ.get("EMBL_FILE_TYPE", "FLATFILE")

    ENA_TOKEN = environ.get("ENA_TOKEN")

    CONSENSUS_MAJOR_FASTA_SUFFIX = environ.get(
        "CONSENSUS_MAJOR_FASTA_SUFFIX", "_consensus_major.fa.gz"
    )
    CONSENSUS_MINOR_FASTA_SUFFIX = environ.get(
        "CONSENSUS_MINOR_FASTA_SUFFIX", "_consensus_minor.fa.gz"
    )
    CONSENSUS_MAJOR_EMBL_SUFFIX = environ.get(
        "CONSENSUS_MAJOR_EMBL_SUFFIX", "_consensus_major.embl.gz"
    )
    CONSENSUS_MINOR_EMBL_SUFFIX = environ.get(
        "CONSENSUS_MINOR_EMBL_SUFFIX", "_consensus_minor.embl.gz"
    )

#---------------------------------
# Utils
# ---------------------------------
def sort_key(x: SampleCount):
    # if the value is None, it is treated as a very small number for sorting purposes
    return x.rpkm_proportions if x.rpkm_proportions is not None else float('-inf')


def build_ser_payload(sample: Sample, sample_counts: Iterable[SampleCount]) -> dict:
    """
    Build the SER payload for a given sample and its SampleCount queryset/list.
    Handles:
    - metadata lookup
    - taxon_id / serotype selection (by highest rpkm_proportions)
    - coverage computation
    - aliases (sample, experiment, run, analysis)
    - CRAM file list
    """
    logger.debug(f"Creating SER payload for {sample}")
    metadata = Metadata.objects.filter(sample=sample).first()
    if metadata is None:
        raise ValueError(f"No metadata found for sample {sample!r}")
    if not sample_counts:
        raise ValueError(f"No sample_counts provided for sample {sample!r}")

    now = dt.datetime.now().strftime("%Y%m%d%H%M%S%f")
    sorted_sample_counts = sorted(sample_counts, key=sort_key, reverse=True)
    top = sorted_sample_counts[0]
    taxon_id = top.substrain.taxon_id
    serotype = top.substrain.serotype

    ent_date = getattr(metadata, "ent_date", None)
    if ent_date is None:
        # fallback: today (local date)
        collection_date = timezone.localdate().strftime("%Y-%m-%d")
        logger.warning(
            f"Metadata.ent_date is NULL for sample {sample.pseudonymized_id}; using today's date ({collection_date}) as 'collection date' for ENA payload."

        )
    else:
        collection_date = ent_date.strftime("%Y-%m-%d")

    geo_location = metadata.prescriber
    coverage = float(top.coverage) + 0.01
    sample_alias = f"revseq_sample_{sample.pseudonymized_id}_{now}"
    experiment_alias = f"revseq_experiment_{sample.pseudonymized_id}_{now}"
    files_qs = File.objects.filter(sample=sample)
    cram_files: list[str] = []
    for file in files_qs:
        if file.type.postfix == ".cram":
            logger.debug(f"Adding CRAM file {file.path} to SER for {sample}.")
            cram_files.append(file.path)

    payload = {
        "template": "default",
        "data": {
            "sample": {
                "title": f"ReVSeq Sample {sample.pseudonymized_id}",
                "alias": sample_alias,
                "host subject id": sample.pseudonymized_id,
                "taxon_id": taxon_id,
                "collection date": collection_date,
                "geographic location (region and locality)": geo_location,
            },
            "experiment": {
                "alias": experiment_alias,
                "sample_alias": sample_alias,
            },
            "run": {
                "alias": f"revseq_run_{sample.pseudonymized_id}_{now}",
                "experiment_alias": experiment_alias,
                "file": cram_files,
            },
            "analysis": {
                "name": f"analysis_{sample.pseudonymized_id}_{now}",
                "coverage": coverage,
                "program": "BCFtools",
            },
        },
        "files": cram_files,
    }

    if serotype:
        payload["data"]["sample"]["serotype"] = serotype
    logger.debug(f"Payload for {sample.pseudonymized_id}: {payload}")
    return payload


def extract_basename(path: str) -> str:
    return os.path.basename(path)

def select_analysis_files_for_sample(
    sample,
    files,
    *,
    consensus_file_suffix: str,
    chromosome_file_name: str,
    embl_file_suffix: str,
    influenza_only: bool,
) -> list:
    """
    Decide which files to upload as analysis input for this sample.
    - If `influenza_only` is True: only EMBL files (matching `embl_file_suffix`).
    - Otherwise: consensus + chromosome list files.
    """
    if influenza_only:
        logger.info(f"We only upload EMBL files for influenza samples like {sample}")

        embl_files = [
            f for f in files
            if extract_basename(f.path).endswith(embl_file_suffix)
        ]

        # TEMP FIX: keep only the "filtered" EMBL if duplicates exist
        if len(embl_files) > 1:
            filtered = [f for f in embl_files if "_filtered" in extract_basename(f.path)]
            if filtered:
                logger.warning(
                    f"Multiple EMBL files found for {sample}. "
                    f"Keeping only filtered one: {filtered[0].path}"
                )
                return filtered

        return embl_files

        # -------------------------
        # Non-influenza: consensus + chromosome
        # -------------------------
    logger.info(f"Uploading consensus and chromosome files for {sample}")
    return [
        f for f in files
        if extract_basename(f.path).endswith(consensus_file_suffix)
           or extract_basename(f.path) == chromosome_file_name
    ]


