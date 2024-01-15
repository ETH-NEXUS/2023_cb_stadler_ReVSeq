from django.core.management.base import BaseCommand


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            '-t', '--type', type=str, choices=['study', 'ser', 'analysis', 'all'],
            help='Type of data to upload: study, ser (sample-experiment-run), analysis, or all'
        )

    def handle(self, *args, **options):
        data_type = options.get('type')

        if data_type == 'study':
            self.upload_study()
        elif data_type == 'ser':
            self.upload_ser()
        elif data_type == 'analysis':
            self.upload_analysis()
        elif data_type == 'all':
            self.upload_study()
            self.upload_ser()
            self.upload_analysis()
        else:
            print('Please specify a data type to upload: study, ser, analysis, or all')

    def upload_study(self):
        pass

    def upload_ser(self):
        pass

    def upload_analysis(self):
        pass
