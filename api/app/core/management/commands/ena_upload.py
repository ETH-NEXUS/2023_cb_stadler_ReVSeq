from django.core.management.base import BaseCommand
import requests
from os import environ
from core.models import Sample, File, Metadata
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


class Command(BaseCommand):
    def __init__(self):
        super().__init__()
        self.job_analysis_ids = []
        self.job_ids = []

    token = environ.get('ENA_TOKEN')
    headers = {'Authorization': f'Token {token}', 'Content-Type': 'application/json'}

    def add_arguments(self, parser):
        parser.add_argument(
            '-t', '--type', type=str, choices=['study', 'ser_and_analysis', 'delete_job_id'],
            help='Type of data to upload: study, or ser_and_analysis to upload the samples which don;t already have a '
                 'job_id'
        )

    def delete_job_id(self):
        samples = Sample.objects.filter(job_id__isnull=False)
        for sample in samples:
            sample.job_id = None
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

    def _release_jobs_loop(self, endpoint):
        try:
            continue_releasing = True
            while continue_releasing:
                jobs = self.handle_http_request(endpoint, method='get')['results']
                submitted_jobs = [j for j in jobs if j['status'] == 'SUBMITTED']
                queued_jobs = [j for j in jobs if j['status'] == 'QUEUED']
                running_jobs = [j for j in jobs if j['status'] == 'RUNNING']
                print(
                    f"There are {len(submitted_jobs)} submitted jobs, {len(queued_jobs)} queued jobs, and {len(running_jobs)} running jobs")
                for job in submitted_jobs:
                    self.release_job(RELEASE_JOB_ENDPOINT, job['id'])
                    print(f"Released job {job['id']}")
                continue_releasing = queued_jobs or running_jobs
                time.sleep(30)
        except Exception as e:
            print(f"An error occurred: {e}")

    def upload_ser(self):
        samples = Sample.objects.filter(job_id__isnull=True, valid=True) # im admin interface add new filed released_to_ena and then filter by this field and sent_to_ena === true
        for sample in samples:
            sample_counts = sample.samplecounts.all()
            if not sample_counts:
                logger.warning(f'No counts for {sample}')
                continue

            payload = self._create_ser_payload(sample, sample_counts)
            print(f"Payload: {payload}")

            response = self.handle_http_request(SER_ENDPOINT, payload, 'post',
                                                message=f'SER for {sample} uploaded successfully')
            files = File.objects.filter(sample=sample)
            # TODO: we need two files here: consensus.fa und chromosome file, both gzipped
            analysis_files = []
            consensus_fa_file = files.filter(type__postfix='.fa').first()
            gzipped_consensus_fa_file = consensus_fa_file
            # if accession id comes, it worked

            job_id = response['id']
            sample.job_id = job_id
            sample.save()
            time.sleep(10)
            self.upload_analysis_job_and_files(job_id, sample, files)
            self.job_ids.append(job_id)

        self._release_jobs_loop(JOBS_ENDPOINT)
        self._release_jobs_loop(ANALYSIS_JOBS_ENDPOINT)

    def _create_ser_payload(self, sample, sample_counts):
        now = dt.datetime.now().strftime('%Y%m%d%H%M%S%f')
        sorted_sample_counts = sorted(sample_counts, key=self.__sort_key, reverse=True)
        taxon_id = sorted_sample_counts[0].substrain.taxon_id
        metadata = Metadata.objects.filter(sample=sample).first()
        collection_date = metadata.ent_date.strftime("%Y-%m-%d")
        geo_location = metadata.prescriber
        coverage = sorted_sample_counts[0].coverage

        sample_alias = f"revseq_sample_{sample.pseudonymized_id}_{now}"
        experiment_alias = f"revseq_experiment_{sample.pseudonymized_id}_{now}"
        _files = []
        files = File.objects.filter(sample=sample)
        for file in files:
            if file.type.postfix == '.cram':
                _files.append(file.path)
                logger.debug(f'Adding {file.path} to SER for {sample}')

        return {
            'template': 'default',
            'data': {
                'sample': {
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

                }

            },
            'files': _files
        }



    def upload_analysis_job_and_files(self, job_id, sample, files):
        payload = {'template': 'default', 'data': {}, 'job': job_id}
        response = self.handle_http_request(ANALYSIS_ENDPOINT, payload, 'post',
                                            message=f'Analysis job for {sample} uploaded successfully')
        analysis_job_id = response['id']
        sample.analysis_job_id = analysis_job_id
        sample.save()
        self.upload_analysis_files(analysis_job_id, files)

        self.enqueue_analysis_job(analysis_job_id)
        self.job_analysis_ids.append(analysis_job_id)

    def upload_analysis_files(self, analysis_job_id, files):
        for file in files:
            # file types: FASTA und CHROMOSOME_LIST
            # TODO: if it is a chromosome file, file type is CHROMOSOME_LIST
            payload = {'job': analysis_job_id, 'file_name': file.path, "file_type": 'FASTA'} #
            self.handle_http_request(ANALYSIS_FILES_ENDPOINT, payload, 'post',
                                     message=f'Analysis file {file} uploaded successfully')

    def enqueue_analysis_job(self, analysis_job_id):
        url = ENQUEUE_ENDPOINT.replace('<job_id>', str(analysis_job_id))
        self.handle_http_request(url, method='get', message=f'Analysis job {analysis_job_id} enqueued')

    def release_job(self, url, job_id):
        _url = url.replace('<job_id>', str(job_id))
        self.handle_http_request(_url, method='get', message=f'Job {job_id} released successfully')

    def handle(self, *args, **options):
        data_type = options.get('type')

        if data_type == 'study':
            self.upload_study()
        elif data_type == 'ser_and_analysis':
            self.upload_ser()

        elif data_type == 'delete_job_id':
            self.delete_job_id()
        else:
            logger.warning('Please specify a data type to upload: study, ser, analysis')
