from django.core.management.base import BaseCommand
from core.models import Sample, File, Plate, FileType
from helpers.color_log import logger
from core.helpers import compute_checksum

file_data = [
    {"plate": "RVSeqPlate3-m2", "sample": "32WNFL", "file": "/data/RVSeqPlate3-m2/sample_m2-32WNFL/m2-32WNFL.embl.gz"},
    {"plate": "RVSeqPlate3-m2", "sample": "375EUk", "file": "/data/RVSeqPlate3-m2/sample_m2-375EUk/m2-375EUk.embl.gz"},
    {"plate": "RVSeqPlate4-m2", "sample": "itJLqU", "file": "/data/RVSeqPlate4-m2/sample_m2-itJLqU/m2-itJLqU.embl.gz"},
    {"plate": "RVSeqPlate4-m2", "sample": "xwkVbK", "file": "/data/RVSeqPlate4-m2/sample_m2-xwkVbK/m2-xwkVbK.embl.gz"},
    {"plate": "RVSeqPlate4-m2", "sample": "RyXauM", "file": "/data/RVSeqPlate4-m2/sample_m2-RyXauM/m2-RyXauM.embl.gz"},
    {"plate": "RVSeqPlate4-m2", "sample": "jtjCjd", "file": "/data/RVSeqPlate4-m2/sample_m2-jtjCjd/m2-jtjCjd.embl.gz"},
    {"plate": "RVSeqPlate4-m2", "sample": "4UF9K6", "file": "/data/RVSeqPlate4-m2/sample_m2-4UF9K6/m2-4UF9K6.embl.gz"},
    {"plate": "RVSeqPlate4-m2", "sample": "CLwScw", "file": "/data/RVSeqPlate4-m2/sample_m2-CLwScw/m2-CLwScw.embl.gz"},
    {"plate": "RVSeqPlate4-m2", "sample": "kwbwf8", "file": "/data/RVSeqPlate4-m2/sample_m2-kwbwf8/m2-kwbwf8.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "m6243p", "file": "/data/RVSeqPlate5-m2/sample_m2-m6243p/m2-m6243p.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "V5iihT", "file": "/data/RVSeqPlate5-m2/sample_m2-V5iihT/m2-V5iihT.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "kBgpZ4", "file": "/data/RVSeqPlate5-m2/sample_m2-kBgpZ4/m2-kBgpZ4.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "x6NRYc", "file": "/data/RVSeqPlate5-m2/sample_m2-x6NRYc/m2-x6NRYc.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "q2oXnd", "file": "/data/RVSeqPlate5-m2/sample_m2-q2oXnd/m2-q2oXnd.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "39kc9U", "file": "/data/RVSeqPlate5-m2/sample_m2-39kc9U/m2-39kc9U.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "xrRPFj", "file": "/data/RVSeqPlate5-m2/sample_m2-xrRPFj/m2-xrRPFj.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "xYnaUy", "file": "/data/RVSeqPlate5-m2/sample_m2-xYnaUy/m2-xYnaUy.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "qtZxkz", "file": "/data/RVSeqPlate5-m2/sample_m2-qtZxkz/m2-qtZxkz.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "VA9Yi5", "file": "/data/RVSeqPlate5-m2/sample_m2-VA9Yi5/m2-VA9Yi5.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "ss4Fh8", "file": "/data/RVSeqPlate5-m2/sample_m2-ss4Fh8/m2-ss4Fh8.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "DwiwwU", "file": "/data/RVSeqPlate5-m2/sample_m2-DwiwwU/m2-DwiwwU.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "9HhvMh", "file": "/data/RVSeqPlate5-m2/sample_m2-9HhvMh/m2-9HhvMh.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "oLXZ9a", "file": "/data/RVSeqPlate5-m2/sample_m2-oLXZ9a/m2-oLXZ9a.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "SQZYHn", "file": "/data/RVSeqPlate5-m2/sample_m2-SQZYHn/m2-SQZYHn.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "SKJgwA", "file": "/data/RVSeqPlate5-m2/sample_m2-SKJgwA/m2-SKJgwA.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "4rvG7f", "file": "/data/RVSeqPlate5-m2/sample_m2-4rvG7f/m2-4rvG7f.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "xHGR9u", "file": "/data/RVSeqPlate5-m2/sample_m2-xHGR9u/m2-xHGR9u.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "eAr6Lo", "file": "/data/RVSeqPlate5-m2/sample_m2-eAr6Lo/m2-eAr6Lo.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "ydSpNQ", "file": "/data/RVSeqPlate5-m2/sample_m2-ydSpNQ/m2-ydSpNQ.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "Sm9GYA", "file": "/data/RVSeqPlate5-m2/sample_m2-Sm9GYA/m2-Sm9GYA.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "opfXjQ", "file": "/data/RVSeqPlate5-m2/sample_m2-opfXjQ/m2-opfXjQ.embl.gz"},
    {"plate": "RVSeqPlate5-m2", "sample": "g33df2", "file": "/data/RVSeqPlate5-m2/sample_m2-g33df2/m2-g33df2.embl.gz"},
    {"plate": "RVSeqPlate6-m2", "sample": "VT3bcF", "file": "/data/RVSeqPlate6-m2/sample_m2-VT3bcF/m2-VT3bcF.embl.gz"},
    {"plate": "RVSeqPlate6-m2", "sample": "mmYEfJ", "file": "/data/RVSeqPlate6-m2/sample_m2-mmYEfJ/m2-mmYEfJ.embl.gz"},
    {"plate": "RVSeqPlate6-m2", "sample": "iDqyV9", "file": "/data/RVSeqPlate6-m2/sample_m2-iDqyV9/m2-iDqyV9.embl.gz"},
    {"plate": "RVSeqPlate6-m2", "sample": "sr4bvM", "file": "/data/RVSeqPlate6-m2/sample_m2-sr4bvM/m2-sr4bvM.embl.gz"},
    {"plate": "RVSeqPlate6-m2", "sample": "LCBdav", "file": "/data/RVSeqPlate6-m2/sample_m2-LCBdav/m2-LCBdav.embl.gz"},
    {"plate": "RVSeqPlate6-m2", "sample": "LosGHn", "file": "/data/RVSeqPlate6-m2/sample_m2-LosGHn/m2-LosGHn.embl.gz"},
    {"plate": "RVSeqPlate6-m2", "sample": "BnbHLr", "file": "/data/RVSeqPlate6-m2/sample_m2-BnbHLr/m2-BnbHLr.embl.gz"},
    {"plate": "RVSeqPlate6-m2", "sample": "uTpHe5", "file": "/data/RVSeqPlate6-m2/sample_m2-uTpHe5/m2-uTpHe5.embl.gz"},
    {"plate": "RVSeqPlate6-m2", "sample": "hzNQRG", "file": "/data/RVSeqPlate6-m2/sample_m2-hzNQRG/m2-hzNQRG.embl.gz"},
    {"plate": "RVSeqPlate6-m2", "sample": "4t7Ncr", "file": "/data/RVSeqPlate6-m2/sample_m2-4t7Ncr/m2-4t7Ncr.embl.gz"},
    {"plate": "RVSeqPlate6-m2", "sample": "K2udH5", "file": "/data/RVSeqPlate6-m2/sample_m2-K2udH5/m2-K2udH5.embl.gz"},
    {"plate": "RVSeqPlate6-m2", "sample": "SpXNmN", "file": "/data/RVSeqPlate6-m2/sample_m2-SpXNmN/m2-SpXNmN.embl.gz"},
    {"plate": "RVSeqPlate6-m2", "sample": "h2ZJd5", "file": "/data/RVSeqPlate6-m2/sample_m2-h2ZJd5/m2-h2ZJd5.embl.gz"},
    {"plate": "RVSeqPlate6-m2", "sample": "juVwzj", "file": "/data/RVSeqPlate6-m2/sample_m2-juVwzj/m2-juVwzj.embl.gz"},
    {"plate": "RVSeqPlate6-m2", "sample": "3H4uLZ", "file": "/data/RVSeqPlate6-m2/sample_m2-3H4uLZ/m2-3H4uLZ.embl.gz"},
    {"plate": "RVSeqPlate6-m2", "sample": "gkcF5C", "file": "/data/RVSeqPlate6-m2/sample_m2-gkcF5C/m2-gkcF5C.embl.gz"},
    {"plate": "RVSeqPlate8-m2", "sample": "gaQkyT", "file": "/data/RVSeqPlate8-m2/sample_m2-gaQkyT/m2-gaQkyT.embl.gz"},
    {"plate": "RVSeqPlate8-m2", "sample": "FXFt6J", "file": "/data/RVSeqPlate8-m2/sample_m2-FXFt6J/m2-FXFt6J.embl.gz"},
    {"plate": "RVSeqPlate8-m2", "sample": "f8rfBe", "file": "/data/RVSeqPlate8-m2/sample_m2-f8rfBe/m2-f8rfBe.embl.gz"},
    {"plate": "RVSeqPlate8-m2", "sample": "uSTYZ2", "file": "/data/RVSeqPlate8-m2/sample_m2-uSTYZ2/m2-uSTYZ2.embl.gz"},
    {"plate": "RVSeqPlate8-m2", "sample": "yhE88h", "file": "/data/RVSeqPlate8-m2/sample_m2-yhE88h/m2-yhE88h.embl.gz"},
    {"plate": "RVSeqPlate8-m2", "sample": "j3ZXpL", "file": "/data/RVSeqPlate8-m2/sample_m2-j3ZXpL/m2-j3ZXpL.embl.gz"},
    {"plate": "RVSeqPlate8-m2", "sample": "2hKeW2", "file": "/data/RVSeqPlate8-m2/sample_m2-2hKeW2/m2-2hKeW2.embl.gz"},
    {"plate": "RVSeqPlate8-m2", "sample": "b99ooj", "file": "/data/RVSeqPlate8-m2/sample_m2-b99ooj/m2-b99ooj.embl.gz"},
    {"plate": "RVSeqPlate8-m2", "sample": "sAyJww", "file": "/data/RVSeqPlate8-m2/sample_m2-sAyJww/m2-sAyJww.embl.gz"},
    {"plate": "RVSeqPlate9-m2", "sample": "sFFUfJ", "file": "/data/RVSeqPlate9-m2/sample_m2-sFFUfJ/m2-sFFUfJ.embl.gz"},
    {"plate": "RVSeqPlate9-m2", "sample": "5cJWyf", "file": "/data/RVSeqPlate9-m2/sample_m2-5cJWyf/m2-5cJWyf.embl.gz"},
    {"plate": "RVSeqPlate9-m2", "sample": "TPgies", "file": "/data/RVSeqPlate9-m2/sample_m2-TPgies/m2-TPgies.embl.gz"},
    {"plate": "RVSeqPlate9-m2", "sample": "6z2wf2", "file": "/data/RVSeqPlate9-m2/sample_m2-6z2wf2/m2-6z2wf2.embl.gz"},
    {"plate": "RVSeqPlate11-m2", "sample": "wVBQbT", "file": "/data/RVSeqPlate11-m2/sample_m2-wVBQbT/m2-wVBQbT.embl.gz"}
]
class Command(BaseCommand):
    help = 'Imports file paths and assigns them to samples from a text file.'

    def handle(self, *args, **kwargs):
        for entry in file_data:
            sample_id = entry["sample"]
            file_path = entry["file"]

            try:
                sample = Sample.objects.get(pseudonymized_id=sample_id)
                file_type, created = FileType.objects.get_or_create(postfix="embl")
                checksum = compute_checksum(file_path)
               # if file with this path dont exist , create
                if File.objects.filter(path=file_path).exists():
                    logger.info(f'File {file_path} already exists for sample {sample_id} on plate {sample.plate}')
                    continue
                else:
                    logger.info(f'File {file_path} does not exist, creating new file instance.')

                    file_instance = File.objects.create(
                        path=file_path,
                        checksum=checksum,
                        type=file_type,
                        sample=sample,
                        plate=sample.plate,
                    )
                    logger.info(f'File {file_path} created with type {file_type.postfix} for sample {sample_id} on plate {sample.plate}')
                    logger.info(f'Added file {file_path} for sample {sample_id}')
            except (Plate.DoesNotExist, Sample.DoesNotExist) as e:
                logger.error(f'Error adding file {file_path} for sample {sample_id}: {e}')