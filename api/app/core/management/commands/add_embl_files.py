from django.core.management.base import BaseCommand
from core.models import Sample, File, Plate, FileType
from helpers.color_log import logger
from core.helpers import compute_checksum

file_data = [
    {"plate": "RVSeqPlate3-m2", "sample": "32WNFL", "file": "/data/revseq/results/gather_results/RVSeqPlate3-m2/sample_m2-32WNFL/32WNFL.embl"},
    {"plate": "RVSeqPlate3-m2", "sample": "375EUk", "file": "/data/revseq/results/gather_results/RVSeqPlate3-m2/sample_m2-375EUk/375EUk.embl"},
    {"plate": "RVSeqPlate4-m2", "sample": "itJLqU", "file": "/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-itJLqU/itJLqU.embl"},
    {"plate": "RVSeqPlate4-m2", "sample": "xwkVbK", "file": "/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-xwkVbK/xwkVbK.embl"},
    {"plate": "RVSeqPlate4-m2", "sample": "RyXauM", "file": "/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-RyXauM/RyXauM.embl"},
    {"plate": "RVSeqPlate4-m2", "sample": "jtjCjd", "file": "/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-jtjCjd/jtjCjd.embl"},
    {"plate": "RVSeqPlate4-m2", "sample": "4UF9K6", "file": "/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-4UF9K6/4UF9K6.embl"},
    {"plate": "RVSeqPlate4-m2", "sample": "CLwScw", "file": "/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-CLwScw/CLwScw.embl"},
    {"plate": "RVSeqPlate4-m2", "sample": "kwbwf8", "file": "/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-kwbwf8/kwbwf8.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "m6243p", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-m6243p/m6243p.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "V5iihT", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-V5iihT/V5iihT.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "kBgpZ4", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-kBgpZ4/kBgpZ4.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "x6NRYc", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-x6NRYc/x6NRYc.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "q2oXnd", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-q2oXnd/q2oXnd.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "39kc9U", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-39kc9U/39kc9U.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "xrRPFj", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-xrRPFj/xrRPFj.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "xYnaUy", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-xYnaUy/xYnaUy.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "qtZxkz", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-qtZxkz/qtZxkz.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "VA9Yi5", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-VA9Yi5/VA9Yi5.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "ss4Fh8", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-ss4Fh8/ss4Fh8.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "DwiwwU", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-DwiwwU/DwiwwU.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "9HhvMh", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-9HhvMh/9HhvMh.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "oLXZ9a", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-oLXZ9a/oLXZ9a.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "SQZYHn", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-SQZYHn/SQZYHn.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "SKJgwA", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-SKJgwA/SKJgwA.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "4rvG7f", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-4rvG7f/4rvG7f.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "xHGR9u", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-xHGR9u/xHGR9u.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "eAr6Lo", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-eAr6Lo/eAr6Lo.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "ydSpNQ", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-ydSpNQ/ydSpNQ.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "Sm9GYA", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-Sm9GYA/Sm9GYA.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "opfXjQ", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-opfXjQ/opfXjQ.embl"},
    {"plate": "RVSeqPlate5-m2", "sample": "g33df2", "file": "/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-g33df2/g33df2.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "VT3bcF", "file": "/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-VT3bcF/VT3bcF.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "mmYEfJ", "file": "/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-mmYEfJ/mmYEfJ.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "iDqyV9", "file": "/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-iDqyV9/iDqyV9.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "sr4bvM", "file": "/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-sr4bvM/sr4bvM.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "LCBdav", "file": "/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-LCBdav/LCBdav.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "LosGHn", "file": "/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-LosGHn/LosGHn.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "BnbHLr", "file": "/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-BnbHLr/BnbHLr.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "uTpHe5", "file": "/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-uTpHe5/uTpHe5.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "hzNQRG", "file": "/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-hzNQRG/hzNQRG.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "4t7Ncr", "file": "/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-4t7Ncr/4t7Ncr.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "K2udH5", "file": "/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-K2udH5/K2udH5.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "SpXNmN", "file": "/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-SpXNmN/SpXNmN.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "h2ZJd5", "file": "/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-h2ZJd5/h2ZJd5.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "juVwzj", "file": "/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-juVwzj/juVwzj.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "3H4uLZ", "file": "/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-3H4uLZ/3H4uLZ.embl"},
    {"plate": "RVSeqPlate6-m2", "sample": "gkcF5C", "file": "/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-gkcF5C/gkcF5C.embl"},
    {"plate": "RVSeqPlate8-m2", "sample": "gaQkyT", "file": "/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-gaQkyT/gaQkyT.embl"},
    {"plate": "RVSeqPlate8-m2", "sample": "FXFt6J", "file": "/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-FXFt6J/FXFt6J.embl"},
    {"plate": "RVSeqPlate8-m2", "sample": "f8rfBe", "file": "/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-f8rfBe/f8rfBe.embl"},
    {"plate": "RVSeqPlate8-m2", "sample": "uSTYZ2", "file": "/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-uSTYZ2/uSTYZ2.embl"},
    {"plate": "RVSeqPlate8-m2", "sample": "yhE88h", "file": "/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-yhE88h/yhE88h.embl"},
    {"plate": "RVSeqPlate8-m2", "sample": "j3ZXpL", "file": "/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-j3ZXpL/j3ZXpL.embl"},
    {"plate": "RVSeqPlate8-m2", "sample": "2hKeW2", "file": "/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-2hKeW2/2hKeW2.embl"},
    {"plate": "RVSeqPlate8-m2", "sample": "b99ooj", "file": "/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-b99ooj/b99ooj.embl"},
    {"plate": "RVSeqPlate8-m2", "sample": "sAyJww", "file": "/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-sAyJww/sAyJww.embl"},
    {"plate": "RVSeqPlate9-m2", "sample": "sFFUfJ", "file": "/data/revseq/results/gather_results/RVSeqPlate9-m2/sample_m2-sFFUfJ/sFFUfJ.embl"},
    {"plate": "RVSeqPlate9-m2", "sample": "5cJWyf", "file": "/data/revseq/results/gather_results/RVSeqPlate9-m2/sample_m2-5cJWyf/5cJWyf.embl"},
    {"plate": "RVSeqPlate9-m2", "sample": "TPgies", "file": "/data/revseq/results/gather_results/RVSeqPlate9-m2/sample_m2-TPgies/TPgies.embl"},
    {"plate": "RVSeqPlate9-m2", "sample": "6z2wf2", "file": "/data/revseq/results/gather_results/RVSeqPlate9-m2/sample_m2-6z2wf2/6z2wf2.embl"},
    {"plate": "RVSeqPlate11-m2", "sample": "wVBQbT", "file": "/data/revseq/results/gather_results/RVSeqPlate11-m2/sample_m2-wVBQbT/wVBQbT.embl"},
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