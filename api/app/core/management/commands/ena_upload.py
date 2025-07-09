import os

from django.core.management.base import BaseCommand
import requests
from os import environ
from core.models import Sample, File, Metadata, SampleCount
import datetime as dt
from helpers.color_log import logger
import time

STUDY_ENDPOINT = 'http://ena:5000/api/jobs/study/'
SER_ENDPOINT = 'http://ena:5000/api/jobs/ser/'
ANALYSIS_ENDPOINT = 'http://ena:5000/api/analysisjobs/'
ANALYSIS_FILES_ENDPOINT = 'http://ena:5000/api/analysisfiles/'
ENQUEUE_ENDPOINT = 'http://ena:5000/api/analysisjobs/<job_id>/enqueue/'
RELEASE_JOB_ENDPOINT = 'http://ena:5000/api/jobs/<job_id>/release/'
RELEASE_ANALYSIS_JOB_ENDPOINT = 'http://ena:5000/api/analysisjobs/<job_id>/release/'
JOBS_ENDPOINT = 'http://ena:5000/api/jobs/'
ANALYSIS_JOBS_ENDPOINT = 'http://ena:5000/api/analysisjobs/'
CONSENSUS_FILE_SUFFIX = '.fa.gz'
CHROMOSOME_FILE_NAME = 'chr_file.txt.gz'
EMBL_FILE_SUFFIX = '.embl.gz'
EMBL_FILE_TYPE = 'FLATFILE'


def extract_basename(path):
    return os.path.basename(path)

class Command(BaseCommand):
    def __init__(self):
        super().__init__()
        self.job_analysis_ids = []
        self.job_ids = []
        self.data_for_analysis_upload = []

    token = environ.get('ENA_TOKEN')
    headers = {'Authorization': f'Token {token}', 'Content-Type': 'application/json'}

    def add_arguments(self, parser):
        parser.add_argument(
            '-t', '--type', type=str, choices=['study', 'ser_and_analysis', 'delete_job_id'],
            help='Type of data to upload: study, or ser_and_analysis to upload the samples which don;t already have a '
                 'job_id'
        )
        # if we have a custom list of samples to upload, we can add it here as a lis []
        parser.add_argument(
            '-s', '--samples', type=str, nargs='*',
            help='List of sample IDs to upload. If not provided, all samples with job_id=None will be uploaded.'
        )

        parser.add_argument(
            '-r', '--task', type=str, choices=['resend_analysis_jobs'],
            help='Task to perform: resend_analysis_jobs to resend analysis jobs for the given samples'
        )
        # argument not to include analysis
        parser.add_argument(
            '-a', '--no_analysis', action='store_true',
            help='If set, analysis jobs will not be uploaded, only SER jobs.'
        )

        # command to upload this one 32WNFL without analysis
        # python manage.py ena_upload --type ser_and_analysis --no_analysis -s 32WNFL
    def resend_analysis_jobs(self, samples):
        if not samples:
            logger.warning('No samples provided for resend_analysis_jobs')
            return
        for sample_id in samples:
            sample = Sample.objects.filter(pseudonymized_id=sample_id).first()
            sample_counts = SampleCount.objects.filter(sample=sample)
            payload = self._create_ser_payload(sample, sample_counts)
            analysis_payload = payload['data']['analysis']
            files = File.objects.filter(sample=sample)
            analysis_files = [file for file in files if extract_basename(file.path).endswith(EMBL_FILE_SUFFIX)]
            job_id = sample.job_id
            self.data_for_analysis_upload.append((job_id, sample, analysis_files, analysis_payload))
        self._upload_analysis_loop()



    def ena_uploadjob_id(self):
        samples = Sample.objects.filter(job_id__isnull=False)
        for sample in samples:
            sample.job_id = None
            sample.analysis_job_id = None
            sample.save()
            logger.info(f'Job id for {sample} deleted')

    def handle_http_request(self, url, payload=None, method='get', message=None):
        try:
            response = None
            if method == 'get':
                response = requests.get(url, headers=self.headers)
            elif method == 'post':
                response = requests.post(url, headers=self.headers, json=payload)
            if response.status_code in [200, 201]:
                if message:
                    logger.info(message)
                return response.json()
            else:
                logger.error(f'Error calling the endpoint {url}: {response.text}')
        except Exception as e:
            logger.error(f'Exception: {e}')
        return None

    def __sort_key(self, x):
        # if the value is None, it is treated as a very small number for sorting purposes
        return x.rpkm_proportions if x.rpkm_proportions is not None else float('-inf')

    def upload_study(self):
        payload = {'template': 'default', 'data': {}}
        self.handle_http_request(STUDY_ENDPOINT, payload, 'post', 'Study uploaded successfully')

    def _release_jobs_loop(self,
                           list_endpoint,
                           release_endpoint_template,
                           job_type_description="jobs"):
        try:
            continue_releasing = True
            while continue_releasing:
                response = self.handle_http_request(list_endpoint, method='get')
                if not response or 'results' not in response:
                    print(f"Failed to retrieve jobs from {list_endpoint}")
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
                    self.release_job(release_endpoint_template, job['id'])


                continue_releasing = bool(queued_jobs or running_jobs)
                time.sleep(30)
        except Exception as e:
            logger.error(f"An error occurred while releasing {job_type_description}: {e}")

    def upload_ser(self, no_analysis=False, given_samples=None):
        if given_samples:
            samples = Sample.objects.filter(pseudonymized_id__in=given_samples)
            print(f'Found {len(samples)} samples for the given pseudonymized IDs: {given_samples}')
        else:
            samples = Sample.objects.filter(job_id__isnull=True, valid=True, upload_to_ena=True)
            if not samples:
                logger.warning(f'No samples found for the given pseudonymized IDs: {given_samples}')
                return
        logger.info(f'Uploading {len(samples)} samples to ENA')
        print(list(samples.values_list('pseudonymized_id', flat=True)))
        for sample in samples:
            logger.info(f' ##################### Uploading {sample} to ENA #######################')
            sample_counts = SampleCount.objects.filter(sample=sample)
            if not sample_counts:
                logger.warning(f'No counts for {sample}')
                continue
            payload = self._create_ser_payload(sample, sample_counts)
            analysis_payload = payload['data']['analysis']
            response = self.handle_http_request(SER_ENDPOINT, payload, 'post',
                                                message=f'SER for {sample} uploaded successfully')
            files = File.objects.filter(sample=sample)
            analysis_files = [file for file in files if  extract_basename(file.path).endswith(EMBL_FILE_SUFFIX)]
            job_id = response['id']
            sample.job_id = job_id
            sample.save()
            time.sleep(20)
            self.data_for_analysis_upload.append((job_id, sample, analysis_files, analysis_payload))
            self.job_ids.append(job_id)
        if not no_analysis:
            logger.info('Uploading analysis jobs and files')
            self._upload_analysis_loop()
        self._release_jobs_loop(
            list_endpoint=JOBS_ENDPOINT,
            release_endpoint_template=RELEASE_JOB_ENDPOINT,
            job_type_description="regular jobs"
        )
        # self._release_jobs_loop(
        #     list_endpoint=ANALYSIS_JOBS_ENDPOINT,
        #     release_endpoint_template=RELEASE_ANALYSIS_JOB_ENDPOINT,
        #     job_type_description="analysis jobs"
        # )

    def _upload_analysis_loop(self):
        logger.debug(f'Uploading analysis jobs for {len(self.data_for_analysis_upload)} samples')
        for job_id, sample, files, analysis_payload in self.data_for_analysis_upload:
            self.upload_analysis_job_and_files(job_id, sample, files, analysis_payload)
            time.sleep(20)

    def _create_ser_payload(self, sample, sample_counts):
        logger.debug(f'Creating SER payload for {sample}')
        now = dt.datetime.now().strftime('%Y%m%d%H%M%S%f')
        sorted_sample_counts = sorted(sample_counts, key=self.__sort_key, reverse=True)
        taxon_id = sorted_sample_counts[0].substrain.taxon_id
        serotype = sorted_sample_counts[0].substrain.serotype
        metadata = Metadata.objects.filter(sample=sample).first()
        collection_date = metadata.ent_date.strftime("%Y-%m-%d")
        geo_location = metadata.prescriber
        coverage = float(sorted_sample_counts[0].coverage) + 0.01
        sample_alias = f"revseq_sample_{sample.pseudonymized_id}_{now}"
        experiment_alias = f"revseq_experiment_{sample.pseudonymized_id}_{now}"
        _files = []
        files = File.objects.filter(sample=sample)
        for file in files:
            if file.type.postfix == '.cram':
                logger.debug(f'Adding CRAM FILE{file.path} to SER for {sample} ___________________________-')
                _files.append(file.path)
        payload = {
            'template': 'default',
            'data': {
                'sample': {
                    'title': f"ReVSeq Sample {sample.pseudonymized_id}",
                    'alias': sample_alias,
                    'host subject id': sample.pseudonymized_id,
                    'taxon_id': taxon_id,
                    'collection date': collection_date,
                    'geographic location (region and locality)': geo_location
                },
                'experiment': {
                    'alias': experiment_alias,
                    'sample_alias': sample_alias
                },
                'run': {
                    'alias': f"revseq_run_{sample.pseudonymized_id}_{now}",
                    'experiment_alias': experiment_alias,
                    'file': _files
                },
                'analysis': {
                    'name': f"analysis_{sample.pseudonymized_id}_{now}",
                    'coverage': coverage,
                    'program': 'BCFtools',
                }

            },
            'files': _files
        }

        if serotype:
            payload['data']['sample']['serotype'] = serotype
        return payload

    def upload_analysis_job_and_files(self, job_id, sample, files, analysis_payload):
        payload = {'template': 'default', 'data': analysis_payload, 'job': job_id}
        logger.info(f'Uploading analysis payload job for {sample}')
        response = self.handle_http_request(ANALYSIS_ENDPOINT, payload, 'post',
                                            message=f'Analysis job for {sample} uploaded successfully')
        if not response:
            logger.error(f'Failed to upload analysis job for {sample}')
            return
        analysis_job_id = response['id']
        sample.analysis_job_id = analysis_job_id
        sample.save()
        self.upload_analysis_files(analysis_job_id, files)
        time.sleep(30)
        self.enqueue_analysis_job(analysis_job_id)
        self.job_analysis_ids.append(analysis_job_id)

    def upload_analysis_files(self, analysis_job_id, files):
        for file in files:
            # file types: FASTA und CHROMOSOME_LIST
            # if it is a chromosome file, file type is CHROMOSOME_LIST
            # file_type = 'FASTA'
            # if extract_basename(file.path) == CHROMOSOME_FILE_NAME:
            #     file_type = 'CHROMOSOME_LIST'
            # payload = {'job': analysis_job_id, 'file_name': file.path, "file_type": file_type} #
            # self.handle_http_request(ANALYSIS_FILES_ENDPOINT, payload, 'post',
            #                          message=f'Analysis file {file} uploaded successfully')

            # for EMBL files, we need to set the file type to FLATFILE
            # we don't need fasta and chromosome files anymore, only EMBL files
            if extract_basename(file.path).endswith(EMBL_FILE_SUFFIX):
                payload = {'job': analysis_job_id, 'file_name': file.path, "file_type": EMBL_FILE_TYPE}
                self.handle_http_request(ANALYSIS_FILES_ENDPOINT, payload, 'post',
                                         message=f'Analysis file {file} uploaded successfully')
            else:
                logger.warning(f'Skipping file {file.path} for analysis upload, not an EMBL file')


    def enqueue_analysis_job(self, analysis_job_id):
        url = ENQUEUE_ENDPOINT.replace('<job_id>', str(analysis_job_id))
        self.handle_http_request(url, method='get', message=f'Analysis job {analysis_job_id} enqueued')

    def release_job(self, url, job_id):
        _url = url.replace('<job_id>', str(job_id))
        self.handle_http_request(_url, method='get', message=f'Job {job_id} released successfully')


    def handle(self, *args, **options):
        task = options.get('task')
        no_analysis = options.get('no_analysis', False)
        samples = options.get('samples', None)
        if task == 'resend_analysis_jobs':
            samples = options.get('samples')
            self.resend_analysis_jobs(samples)
            return
        data_type = options.get('type')
        if data_type == 'study':
            self.upload_study()
        elif data_type == 'ser_and_analysis':
            self.upload_ser(no_analysis, samples)
        elif data_type == 'delete_job_id':
            self.delete_job_id()

        else:
            logger.warning('Please specify a data type to upload: study, ser, analysis')
