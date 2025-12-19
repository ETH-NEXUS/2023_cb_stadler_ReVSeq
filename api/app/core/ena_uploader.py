from typing import List
from helpers.color_log import logger
from .utils import handle_http_request, Defaults, build_ser_payload, select_analysis_files_for_sample, extract_basename
from core.models import Sample, SampleCount, File, Metadata
from copy import deepcopy
import time

ANALYSIS_KIND_STANDARD = "standard"
ANALYSIS_KIND_COINF_MAJOR = "coinf_major"
ANALYSIS_KIND_COINF_MINOR = "coinf_minor"

class ENAUploader:
    def __init__(
            self,
            ena_token: str | None = Defaults.ENA_TOKEN,
            *,  # this means all following parameters must be named
            study_endpoint: str = Defaults.STUDY_ENDPOINT,
            ser_endpoint: str = Defaults.SER_ENDPOINT,
            analysis_endpoint: str = Defaults.ANALYSIS_ENDPOINT,
            analysis_files_endpoint: str = Defaults.ANALYSIS_FILES_ENDPOINT,
            enqueue_endpoint: str = Defaults.ENQUEUE_ENDPOINT,
            release_job_endpoint: str = Defaults.RELEASE_JOB_ENDPOINT,
            release_analysis_job_endpoint: str = Defaults.RELEASE_ANALYSIS_JOB_ENDPOINT,
            jobs_endpoint: str = Defaults.JOBS_ENDPOINT,
            analysis_jobs_endpoint: str = Defaults.ANALYSIS_JOBS_ENDPOINT,
            consensus_file_suffix: str = Defaults.CONSENSUS_FILE_SUFFIX,
            chromosome_file_name: str = Defaults.CHROMOSOME_FILE_NAME,
            embl_file_suffix: str = Defaults.EMBL_FILE_SUFFIX,
            embl_file_type: str = Defaults.EMBL_FILE_TYPE,
            consensus_major_fasta_suffix: str = Defaults.CONSENSUS_MAJOR_FASTA_SUFFIX,
            consensus_minor_fasta_suffix: str = Defaults.CONSENSUS_MINOR_FASTA_SUFFIX,
            consensus_major_embl_suffix: str = Defaults.CONSENSUS_MAJOR_EMBL_SUFFIX,
            consensus_minor_embl_suffix: str = Defaults.CONSENSUS_MINOR_EMBL_SUFFIX,
            chromosome_major_file_suffix: str = Defaults.CHROMOSOME_MAJOR_FILE_SUFFIX,
            chromosome_minor_file_suffix: str = Defaults.CHROMOSOME_MINOR_FILE_SUFFIX,
            test_run: bool = True


    ):
        if ena_token is None:
            raise ValueError("ENA token is not configured (ENA_TOKEN is missing).")

        self.token = ena_token
        self.headers = {
            'Authorization': f'Token {self.token}',
            'Content-Type': 'application/json',
        }
        # endpoints
        self.study_endpoint = study_endpoint
        self.ser_endpoint = ser_endpoint
        self.analysis_endpoint = analysis_endpoint
        self.analysis_files_endpoint = analysis_files_endpoint
        self.enqueue_endpoint = enqueue_endpoint
        self.release_job_endpoint = release_job_endpoint
        self.release_analysis_job_endpoint = release_analysis_job_endpoint
        self.jobs_endpoint = jobs_endpoint
        self.analysis_jobs_endpoint = analysis_jobs_endpoint
        # file naming / types
        self.consensus_file_suffix = consensus_file_suffix
        self.chromosome_major_file_suffix = chromosome_major_file_suffix
        self.chromosome_minor_file_suffix = chromosome_minor_file_suffix
        self.chromosome_file_name = chromosome_file_name
        self.embl_file_suffix = embl_file_suffix
        self.embl_file_type = embl_file_type
        self.chromosome_file_type = 'CHROMOSOME_LIST'
        self.fasta_file_type = 'FASTA'
        self.consensus_major_fasta_suffix = consensus_major_fasta_suffix
        self.consensus_minor_fasta_suffix = consensus_minor_fasta_suffix
        self.consensus_major_embl_suffix = consensus_major_embl_suffix
        self.consensus_minor_embl_suffix = consensus_minor_embl_suffix
        # behaviour flags
        self.test_run = test_run
        # internal state
        self.job_ids: list[int] = []
        self.job_analysis_ids: list[int] = []
        self.data_for_analysis_upload: list[tuple] = []


    def upload_study(self):
        payload = {'template': 'default', 'data': {}}
        handle_http_request(self.study_endpoint, payload, 'post', 'Study uploaded successfully', self.headers)

    def upload_ser(
        self,
        *,
        pseudonymized_ids: List[str],
        no_analysis: bool = False,
        influenza_only: bool = False,
        coinfections: bool = False,
    ) -> None:
        """
        Upload SER (and optionally analysis) for the given samples.
        """
        if not pseudonymized_ids:
            logger.warning('No pseudonymized IDs provided for upload_ser')
            return
        samples = Sample.objects.filter(pseudonymized_id__in=pseudonymized_ids)
        logger.info(f'Uploading {len(samples)} samples to ENA')
        print(list(samples.values_list('pseudonymized_id', flat=True)))
        if not samples.exists():
            logger.warning(f'No samples found for the given pseudonymized IDs: {pseudonymized_ids}')
            return
        for sample in samples:
            self.upload_single_sample(
                sample,
                influenza_only=influenza_only,
                coinfections=coinfections,
            )

        if not no_analysis:
            logger.info('**************  Uploading analysis jobs and files ***************')
            self.upload_analysis_loop()

        self.release_jobs_loop(
            list_endpoint=self.jobs_endpoint,
            release_endpoint_url_template=self.release_job_endpoint,
            job_type_description="regular jobs"
        )

    # ------------------------- Single sample upload ------------------------- #

    def upload_single_sample(
            self,
            sample: Sample,
            *,
            influenza_only: bool,
            coinfections: bool = False,
    ) -> None:
        logger.info(
            f"--------------------------- Uploading {sample} to ENA -----------------------------"
        )
        sample_counts = SampleCount.objects.filter(sample=sample)
        if not sample_counts:
            logger.warning("No counts for %s", sample)
            return

        payload = build_ser_payload(sample, sample_counts)
        analysis_payload = payload["data"]["analysis"]

        response = handle_http_request(
            self.ser_endpoint,
            payload,
            "post",
            message=f"SER for {sample} uploaded successfully",
            headers=self.headers,
        )
        if not response:
            logger.error("Failed to upload SER for %s", sample)
            return
        job_id = response["id"]
        if self.test_run:
            logger.info(
                "Test run enabled - saving test_job_id %s to sample %s",
                job_id,
                sample,
            )
            sample.test_job_id = job_id
        else:
            sample.job_id = job_id
        sample.save()

        files = File.objects.filter(sample=sample)

        if coinfections:
            logger.info(f'Preparing co-infection analysis upload for sample {sample} (job_id={job_id})')
            self._prepare_coinfection_analysis_upload(
                job_id=job_id,
                sample=sample,
                files=files,
                analysis_payload=analysis_payload,
            )
        else:
            self._prepare_standard_analysis_upload(
                job_id=job_id,
                sample=sample,
                files=files,
                analysis_payload=analysis_payload,
                influenza_only=influenza_only,
            )
        self.job_ids.append(job_id)

        # ------------------------- Analysis upload loop and related methods ------------------------- #
    def _prepare_standard_analysis_upload(
            self,
            *,
            job_id,
            sample,
            files,
            analysis_payload,
            influenza_only: bool,
    ) -> None:
        analysis_files = select_analysis_files_for_sample(
            sample,
            files,
            consensus_file_suffix=self.consensus_file_suffix,
            chromosome_file_name=self.chromosome_file_name,
            embl_file_suffix=self.embl_file_suffix,
            influenza_only=influenza_only,
        )
        if not analysis_files:
            logger.warning(
                f'!!!!!!!!!!!!!!!! No analysis files selected for sample {sample}; skipping analysis upload. !!!!!!!!!!!!!!!'
            )
            return

        self.data_for_analysis_upload.append(
            (job_id, sample, analysis_files, analysis_payload, ANALYSIS_KIND_STANDARD)
        )

    def _prepare_coinfection_analysis_upload(
            self,
            *,
            job_id: int,
            sample: Sample,
            files,
            analysis_payload: dict,
    ) -> None:
        """
        Co-infections:
          - one SER job_id
          - up to two analysis jobs: MAJOR + MINOR

        For each side independently:
          - consensus: prefer EMBL if present, else FASTA
          - chromosome list: attach matching chr file if present
          - if a side is missing consensus, skip that side (log)
          - chr file is optional: if missing, still upload consensus (log)
        """

        def bn(f: File) -> str:
            return extract_basename(f.path)

        # ---- collect candidates per side ----
        major_embl = major_fa = None
        minor_embl = minor_fa = None
        major_chr = None
        minor_chr = None

        for f in files:
            name = bn(f)

            # MAJOR consensus
            if name.endswith(self.consensus_major_embl_suffix):
                major_embl = f
                continue
            if name.endswith(self.consensus_major_fasta_suffix):
                major_fa = f
                continue

            # MINOR consensus
            if name.endswith(self.consensus_minor_embl_suffix):
                minor_embl = f
                continue
            if name.endswith(self.consensus_minor_fasta_suffix):
                minor_fa = f
                continue

            # chr files (suffix-based)
            if name.endswith(self.chromosome_major_file_suffix):
                major_chr = f
                continue
            if name.endswith(self.chromosome_minor_file_suffix):
                minor_chr = f
                continue

        # choose best consensus per side
        major_consensus = major_embl or major_fa
        minor_consensus = minor_embl or minor_fa

        if major_consensus is None and minor_consensus is None:
            logger.warning(
                f"Co-infection {sample} (job_id={job_id}): no major/minor consensus found; skipping analysis upload."
            )
            return

        base_name = analysis_payload.get("name", f"analysis_{sample.pseudonymized_id}")

        # ---- MAJOR ----
        if major_consensus is not None:
            major_payload = deepcopy(analysis_payload)
            major_payload["name"] = f"{base_name}_major"

            major_files: list[File] = [major_consensus]
            if major_chr is not None:
                major_files.append(major_chr)
                logger.info(
                    f"Co-infection {sample} (job_id={job_id}): MAJOR will upload consensus={major_consensus.path} + chr={major_chr.path}"
                )
            else:
                logger.warning(
                    f"Co-infection {sample} (job_id={job_id}): MAJOR chr file missing (suffix={self.chromosome_major_file_suffix}); uploading consensus only ({major_consensus.path})"
                )

            self.data_for_analysis_upload.append(
                (job_id, sample, major_files, major_payload, ANALYSIS_KIND_COINF_MAJOR)
            )
        else:
            logger.warning(
                f"Co-infection {sample} (job_id={job_id}): MAJOR consensus missing; skipping MAJOR analysis upload."
            )

        # ---- MINOR ----
        if minor_consensus is not None:
            minor_payload = deepcopy(analysis_payload)
            minor_payload["name"] = f"{base_name}_minor"

            minor_files: list[File] = [minor_consensus]
            if minor_chr is not None:
                minor_files.append(minor_chr)
                logger.info(
                    f"Co-infection {sample} (job_id={job_id}): MINOR will upload consensus={minor_consensus.path} + chr={minor_chr.path}"
                )
            else:
                logger.warning(
                    f"Co-infection {sample} (job_id={job_id}): MINOR chr file missing (suffix={self.chromosome_minor_file_suffix}); uploading consensus only ({minor_consensus.path})"
                )

            self.data_for_analysis_upload.append(
                (job_id, sample, minor_files, minor_payload, ANALYSIS_KIND_COINF_MINOR)
            )
        else:
            logger.warning(
                f"Co-infection {sample} (job_id={job_id}): MINOR consensus missing; skipping MINOR analysis upload."
            )

    def upload_analysis_loop(self) -> None:
        logger.debug(f'Uploading analysis jobs for {len(self.data_for_analysis_upload)} samples')
        for job_id, sample, analysis_files, analysis_payload, kind in self.data_for_analysis_upload:
            self._upload_analysis_job_and_files(job_id, sample, analysis_files, analysis_payload, kind)
            time.sleep(20)

    def _upload_analysis_job_and_files(self, job_id, sample, analysis_files, analysis_payload, kind: str) -> None:
        payload = {'template': 'default', 'data': analysis_payload, 'job': job_id}
        logger.info(f'Uploading analysis payload job for {sample}')
        logger.debug(f'Analysis payload: {payload}')
        logger.info(f'Files for analysis upload: {[file.path for file in analysis_files]}')
        response = handle_http_request(self.analysis_endpoint, payload, 'post',
                                       message=f'Analysis job for {sample} uploaded successfully',
                                        headers=self.headers
                                       )
        if not response:
            logger.error(f'Failed to upload analysis job for {sample}')
            return

        analysis_job_id = response['id']

        if self.test_run:
            # --- TEST RUN FIELDS ---
            if kind == ANALYSIS_KIND_STANDARD:
                sample.test_analysis_job_id = analysis_job_id
            elif kind == ANALYSIS_KIND_COINF_MAJOR:
                sample.test_coinfections_major_analysis_job_id = analysis_job_id
            elif kind == ANALYSIS_KIND_COINF_MINOR:
                sample.test_coinfections_minor_analysis_job_id = analysis_job_id
            else:
                logger.warning(
                    "Unknown analysis kind=%s for sample %s; not saving TEST analysis job id.",
                    kind,
                    sample,
                )
        else:
            # --- PRODUCTION FIELDS ---
            if kind == ANALYSIS_KIND_STANDARD:
                sample.analysis_job_id = analysis_job_id
            elif kind == ANALYSIS_KIND_COINF_MAJOR:
                sample.coinfections_major_analysis_job_id = analysis_job_id
            elif kind == ANALYSIS_KIND_COINF_MINOR:
                sample.coinfections_minor_analysis_job_id = analysis_job_id
            else:
                logger.warning(
                    "Unknown analysis kind=%s for sample %s; not saving analysis job id.",
                    kind,
                    sample,
                )

        sample.save()

        self._upload_analysis_files(analysis_job_id, analysis_files)
        time.sleep(30) # we have to wait a bit before enqueuing
        self._enqueue_analysis_job(analysis_job_id)
        self.job_analysis_ids.append(analysis_job_id)

    def _upload_analysis_files(self, analysis_job_id, analysis_files) -> None:
        for file in analysis_files:
            base = extract_basename(file.path)

            file_type = self.fasta_file_type  # default

            # chromosome lists: classic + coinfections major/minor
            if (
                    base == self.chromosome_file_name
                    or base.endswith(self.chromosome_major_file_suffix)
                    or base.endswith(self.chromosome_minor_file_suffix)
            ):
                file_type = self.chromosome_file_type

            # EMBL flatfiles
            if base.endswith(self.embl_file_suffix):
                file_type = self.embl_file_type

            payload = {"job": analysis_job_id, "file_name": file.path, "file_type": file_type}
            handle_http_request(
                self.analysis_files_endpoint,
                payload,
                "post",
                message=f"Analysis file {file} uploaded successfully",
                headers=self.headers,
            )

    def _enqueue_analysis_job(self, analysis_job_id) -> None:
        url = self.enqueue_endpoint.replace('<job_id>', str(analysis_job_id))
        handle_http_request(
            url,
            method='get',
            message=f'Analysis job {analysis_job_id} enqueued',
            headers=self.headers,
        )



    #  ------------------------- Releasing jobs loop ------------------------- #

    def release_jobs_loop(self,
                           list_endpoint,
                           release_endpoint_url_template,
                           job_type_description="jobs"):
        try:
            continue_releasing = True
            while continue_releasing:
                response = handle_http_request(list_endpoint, method='get', headers=self.headers)
                if not response or 'results' not in response:
                    logger.error(f"Failed to retrieve jobs from {list_endpoint}")
                    break
                jobs = response['results']
                submitted_jobs = [j for j in jobs if j['status'] == 'SUBMITTED']
                queued_jobs = [j for j in jobs if j['status'] == 'QUEUED']
                running_jobs = [j for j in jobs if j['status'] == 'RUNNING']
                error_jobs = [j for j in jobs if j['status'] == 'ERROR']
                logger.info(f"Releasing {job_type_description}: There are {len(submitted_jobs)} submitted, "
                            f"{len(queued_jobs)} queued, and {len(running_jobs)} running jobs.")
                logger.debug(f"Releasing {job_type_description}: There are {len(error_jobs)} error jobs.")

                for job in submitted_jobs:
                    self._release_job(release_endpoint_url_template, job['id'])
                continue_releasing = bool(queued_jobs or running_jobs)
                time.sleep(30)
        except Exception as e:
            logger.error(f"An error occurred while releasing {job_type_description}: {e}")

    def _release_job(self, release_endpoint_url_template, job_id):
        _url = release_endpoint_url_template.replace('<job_id>', str(job_id))
        handle_http_request(_url, method='get', message=f'Job {job_id} released successfully', headers=self.headers)

    # ---------------------------- resend analysis jobs --------------------------- #

    def resend_analysis_jobs(
            self,
            *,
            pseudonymized_ids: List[str],
            influenza_only: bool = False,
    ) -> None:
        """
        Re-send analysis jobs for samples that already have a SER job_id.
        This:
        - rebuilds the analysis payload from current DB state
        - re-selects analysis input files
        - creates new analysis jobs and enqueues them
        """
        if not pseudonymized_ids:
            logger.warning("No pseudonymized IDs provided for resend_analysis_jobs")
            return
        samples = Sample.objects.filter(pseudonymized_id__in=pseudonymized_ids)
        if not samples.exists():
            logger.warning(
                "No samples found for the given pseudonymized IDs: %s",
                pseudonymized_ids,
            )
            return
        logger.info(f"Resending analysis jobs for {samples.count()} samples",)
        logger.debug(f"Resolved samples for resend_analysis_jobs: {list(samples.values_list('pseudonymized_id', flat=True))}")
        # reset per-run state
        self.job_analysis_ids.clear()
        self.data_for_analysis_upload.clear()
        for sample in samples:
            if self.test_run and getattr(sample, "test_job_id", None) is not None:
                job_id = sample.test_job_id
            else:
                job_id = sample.job_id
            if job_id is None:
                logger.warning(
                    f"Sample {sample} has no existing job_id (or test_job_id); "
                    "cannot resend analysis.",
                )
                continue
            logger.info(
                "Preparing analysis re-upload for sample %s (job_id=%s)",
                sample,
                job_id,
            )
            sample_counts = SampleCount.objects.filter(sample=sample)
            if not sample_counts:
                logger.warning("No counts for %s", sample)
                continue
            # build SER payload, extract analysis part
            payload = build_ser_payload(sample, sample_counts)
            analysis_payload = payload["data"]["analysis"]
            files = File.objects.filter(sample=sample)
            analysis_files = select_analysis_files_for_sample(
                sample,
                files,
                consensus_file_suffix=self.consensus_file_suffix,
                chromosome_file_name=self.chromosome_file_name,
                embl_file_suffix=self.embl_file_suffix,
                influenza_only=influenza_only,
            )
            if not analysis_files:
                logger.warning(
                    f"No analysis files selected for sample {sample}; skipping.",
                    sample,
                )
                continue
            # add to the same queue structure used by upload_ser
            self.data_for_analysis_upload.append(
                (job_id, sample, analysis_files, analysis_payload, ANALYSIS_KIND_STANDARD)
            )
        if not self.data_for_analysis_upload:
            logger.warning("No analysis jobs prepared for resend; nothing to do.")
            return
        logger.info(
            f"**************  Resending analysis jobs and files for {len(self.data_for_analysis_upload)} samples ***************"

        )
        self.upload_analysis_loop()

        # ----------------------- co-ninfections ---------------------------- #



