from django.core.management.base import BaseCommand
from core.models import Sample, File, Plate, FileType
from helpers.color_log import logger
from core.helpers import compute_checksum

file_data =[
    {"plate": "RVSeqPlate3-m2", "sample": "m2-32WNFL", "file": "/data/RVSeqPlate3-m2/sample_m2-32WNFL/m2-32WNFL.embl"},
    {"plate": "RVSeqPlate3-m2", "sample": "m2-375EUk", "file": "/data/RVSeqPlate3-m2/sample_m2-375EUk/m2-375EUk.embl"},
    {"plate": "RVSeqPlate4-m2", "sample": "m2-itJLqU", "file": "/data/RVSeqPlate4-m2/sample_m2-itJLqU/m2-itJLqU.embl"},
    {"plate": "RVSeqPlate4-m2", "sample": "m2-xwkVbK", "file": "/data/RVSeqPlate4-m2/sample_m2-xwkVbK/m2-xwkVbK.embl"},
    {"plate": "RVSeqPlate4-m2", "sample": "m2-RyXauM", "file": "/data/RVSeqPlate4-m2/sample_m2-RyXauM/m2-RyXauM.embl"},
    {"plate": "RVSeqPlate4-m2", "sample": "m2-jtjCjd", "file": "/data/RVSeqPlate4-m2/sample_m2-jtjCjd/m2-jtjCjd.embl"},
    {"plate": "RVSeqPlate4-m2", "sample": "m2-4UF9K6", "file": "/data/RVSeqPlate4-m2/sample_m2-4UF9K6/m2-4UF9K6.embl"},
    {"plate": "RVSeqPlate4-m2", "sample": "m2-CLwScw", "file": "/data/RVSeqPlate4-m2/sample_m2-CLwScw/m2-CLwScw.embl"},
    {"plate": "RVSeqPlate4-m2", "sample": "m2-kwbwf8", "file": "/data/RVSeqPlate4-m2/sample_m2-kwbwf8/m2-kwbwf8.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-m6243p", "file": "/data/RVSeqPlate5-m2/sample_m2-m6243p/m2-m6243p.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-V5iihT", "file": "/data/RVSeqPlate5-m2/sample_m2-V5iihT/m2-V5iihT.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-kBgpZ4", "file": "/data/RVSeqPlate5-m2/sample_m2-kBgpZ4/m2-kBgpZ4.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-x6NRYc", "file": "/data/RVSeqPlate5-m2/sample_m2-x6NRYc/m2-x6NRYc.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-q2oXnd", "file": "/data/RVSeqPlate5-m2/sample_m2-q2oXnd/m2-q2oXnd.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-39kc9U", "file": "/data/RVSeqPlate5-m2/sample_m2-39kc9U/m2-39kc9U.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-xrRPFj", "file": "/data/RVSeqPlate5-m2/sample_m2-xrRPFj/m2-xrRPFj.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-xYnaUy", "file": "/data/RVSeqPlate5-m2/sample_m2-xYnaUy/m2-xYnaUy.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-qtZxkz", "file": "/data/RVSeqPlate5-m2/sample_m2-qtZxkz/m2-qtZxkz.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-VA9Yi5", "file": "/data/RVSeqPlate5-m2/sample_m2-VA9Yi5/m2-VA9Yi5.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-ss4Fh8", "file": "/data/RVSeqPlate5-m2/sample_m2-ss4Fh8/m2-ss4Fh8.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-DwiwwU", "file": "/data/RVSeqPlate5-m2/sample_m2-DwiwwU/m2-DwiwwU.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-9HhvMh", "file": "/data/RVSeqPlate5-m2/sample_m2-9HhvMh/m2-9HhvMh.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-oLXZ9a", "file": "/data/RVSeqPlate5-m2/sample_m2-oLXZ9a/m2-oLXZ9a.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-SQZYHn", "file": "/data/RVSeqPlate5-m2/sample_m2-SQZYHn/m2-SQZYHn.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-SKJgwA", "file": "/data/RVSeqPlate5-m2/sample_m2-SKJgwA/m2-SKJgwA.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-4rvG7f", "file": "/data/RVSeqPlate5-m2/sample_m2-4rvG7f/m2-4rvG7f.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-xHGR9u", "file": "/data/RVSeqPlate5-m2/sample_m2-xHGR9u/m2-xHGR9u.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-eAr6Lo", "file": "/data/RVSeqPlate5-m2/sample_m2-eAr6Lo/m2-eAr6Lo.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-ydSpNQ", "file": "/data/RVSeqPlate5-m2/sample_m2-ydSpNQ/m2-ydSpNQ.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-Sm9GYA", "file": "/data/RVSeqPlate5-m2/sample_m2-Sm9GYA/m2-Sm9GYA.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-opfXjQ", "file": "/data/RVSeqPlate5-m2/sample_m2-opfXjQ/m2-opfXjQ.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m2-g33df2", "file": "/data/RVSeqPlate5-m2/sample_m2-g33df2/m2-g33df2.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "m2-VT3bcF", "file": "/data/RVSeqPlate6-m2/sample_m2-VT3bcF/m2-VT3bcF.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "m2-mmYEfJ", "file": "/data/RVSeqPlate6-m2/sample_m2-mmYEfJ/m2-mmYEfJ.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "m2-iDqyV9", "file": "/data/RVSeqPlate6-m2/sample_m2-iDqyV9/m2-iDqyV9.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "m2-sr4bvM", "file": "/data/RVSeqPlate6-m2/sample_m2-sr4bvM/m2-sr4bvM.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "m2-LCBdav", "file": "/data/RVSeqPlate6-m2/sample_m2-LCBdav/m2-LCBdav.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "m2-LosGHn", "file": "/data/RVSeqPlate6-m2/sample_m2-LosGHn/m2-LosGHn.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "m2-BnbHLr", "file": "/data/RVSeqPlate6-m2/sample_m2-BnbHLr/m2-BnbHLr.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "m2-uTpHe5", "file": "/data/RVSeqPlate6-m2/sample_m2-uTpHe5/m2-uTpHe5.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "m2-hzNQRG", "file": "/data/RVSeqPlate6-m2/sample_m2-hzNQRG/m2-hzNQRG.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "m2-4t7Ncr", "file": "/data/RVSeqPlate6-m2/sample_m2-4t7Ncr/m2-4t7Ncr.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "m2-K2udH5", "file": "/data/RVSeqPlate6-m2/sample_m2-K2udH5/m2-K2udH5.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "m2-SpXNmN", "file": "/data/RVSeqPlate6-m2/sample_m2-SpXNmN/m2-SpXNmN.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "m2-h2ZJd5", "file": "/data/RVSeqPlate6-m2/sample_m2-h2ZJd5/m2-h2ZJd5.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "m2-juVwzj", "file": "/data/RVSeqPlate6-m2/sample_m2-juVwzj/m2-juVwzj.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "m2-3H4uLZ", "file": "/data/RVSeqPlate6-m2/sample_m2-3H4uLZ/m2-3H4uLZ.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "m2-gkcF5C", "file": "/data/RVSeqPlate6-m2/sample_m2-gkcF5C/m2-gkcF5C.embl"},
    {"plate": "RVSeqPlate8-m2", "sample": "m2-gaQkyT", "file": "/data/RVSeqPlate8-m2/sample_m2-gaQkyT/m2-gaQkyT.embl"},
    {"plate": "RVSeqPlate8-m2", "sample": "m2-FXFt6J", "file": "/data/RVSeqPlate8-m2/sample_m2-FXFt6J/m2-FXFt6J.embl"},
    {"plate": "RVSeqPlate8-m2", "sample": "m2-f8rfBe", "file": "/data/RVSeqPlate8-m2/sample_m2-f8rfBe/m2-f8rfBe.embl"},
    {"plate": "RVSeqPlate8-m2", "sample": "m2-uSTYZ2", "file": "/data/RVSeqPlate8-m2/sample_m2-uSTYZ2/m2-uSTYZ2.embl"},
    {"plate": "RVSeqPlate8-m2", "sample": "m2-yhE88h", "file": "/data/RVSeqPlate8-m2/sample_m2-yhE88h/m2-yhE88h.embl"},
    {"plate": "RVSeqPlate8-m2", "sample": "m2-j3ZXpL", "file": "/data/RVSeqPlate8-m2/sample_m2-j3ZXpL/m2-j3ZXpL.embl"},
    {"plate": "RVSeqPlate8-m2", "sample": "m2-2hKeW2", "file": "/data/RVSeqPlate8-m2/sample_m2-2hKeW2/m2-2hKeW2.embl"},
    {"plate": "RVSeqPlate8-m2", "sample": "m2-b99ooj", "file": "/data/RVSeqPlate8-m2/sample_m2-b99ooj/m2-b99ooj.embl"},
    {"plate": "RVSeqPlate8-m2", "sample": "m2-sAyJww", "file": "/data/RVSeqPlate8-m2/sample_m2-sAyJww/m2-sAyJww.embl"},
    {"plate": "RVSeqPlate9-m2", "sample": "m2-sFFUfJ", "file": "/data/RVSeqPlate9-m2/sample_m2-sFFUfJ/m2-sFFUfJ.embl"},
    {"plate": "RVSeqPlate9-m2", "sample": "m2-5cJWyf", "file": "/data/RVSeqPlate9-m2/sample_m2-5cJWyf/m2-5cJWyf.embl"},
    {"plate": "RVSeqPlate9-m2", "sample": "m2-TPgies", "file": "/data/RVSeqPlate9-m2/sample_m2-TPgies/m2-TPgies.embl"},
    {"plate": "RVSeqPlate9-m2", "sample": "m2-6z2wf2", "file": "/data/RVSeqPlate9-m2/sample_m2-6z2wf2/m2-6z2wf2.embl"},
    {"plate": "RVSeqPlate11-m2", "sample": "m2-wVBQbT", "file": "/data/RVSeqPlate11-m2/sample_m2-wVBQbT/m2-wVBQbT.embl"}
]
class Command(BaseCommand):
    help = 'Imports file paths and assigns them to samples from a text file.'

    def handle(self, *args, **kwargs):
        for entry in file_data:
            plate_barcode = entry["plate"]
            sample_id = entry["sample"]
            file_path = entry["file"]

            try:
                plate = Plate.objects.get(barcode=plate_barcode)
                sample = Sample.objects.get(pseudonymized_id=sample_id)
                file_type, created = FileType.objects.get_or_create(postfix="embl")
                checksum = compute_checksum(file_path)
               # if file with this path dont exist , create
                if File.objects.filter(path=file_path).exists():
                    logger.info(f'File {file_path} already exists for sample {sample_id} on plate {plate_barcode}')
                    continue
                else:
                    logger.info(f'File {file_path} does not exist, creating new file instance.')

                    file_instance = File.objects.create(
                        path=file_path,
                        checksum=checksum,
                        type=file_type,
                        sample=sample,
                        plate=plate
                    )
                    logger.info(f'File {file_path} created with type {file_type.postfix} for sample {sample_id} on plate {plate_barcode}')
                    logger.info(f'Added file {file_path} for sample {sample_id} on plate {plate_barcode}')
            except (Plate.DoesNotExist, Sample.DoesNotExist) as e:
                logger.error(f'Error adding file {file_path} for sample {sample_id} on plate {plate_barcode}: {e}')