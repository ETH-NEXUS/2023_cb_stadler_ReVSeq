from django.core.management.base import BaseCommand
from django.db import transaction
from core.models import File, Sample, FileType
import hashlib
from pathlib import Path

data = [{'name': 'm2-3qDrEd', 'embl': '/data/revseq/results/gather_results/010623RVS-m2/sample_m2-3qDrEd/m2-3qDrEd.embl.gz'}, {'name': 'm2-4SnyCD', 'embl': '/data/revseq/results/gather_results/010623RVS-m2/sample_m2-4SnyCD/m2-4SnyCD.embl.gz'}, {'name': 'm2-32WNFL', 'embl': '/data/revseq/results/gather_results/RVSeqPlate3-m2/sample_m2-32WNFL/m2-32WNFL.embl.gz'}, {'name': 'm2-BdtP8F', 'embl': '/data/revseq/results/gather_results/RVSeqPlate3-m2/sample_m2-BdtP8F/m2-BdtP8F.embl.gz'}, {'name': 'm2-3URpmc', 'embl': '/data/revseq/results/gather_results/RVSeqPlate3-m2/sample_m2-3URpmc/m2-3URpmc.embl.gz'}, {'name': 'm2-3MgVni', 'embl': '/data/revseq/results/gather_results/RVSeqPlate3-m2/sample_m2-3MgVni/m2-3MgVni.embl.gz'}, {'name': 'm2-4Rt9sw', 'embl': '/data/revseq/results/gather_results/RVSeqPlate3-m2/sample_m2-4Rt9sw/m2-4Rt9sw.embl.gz'}, {'name': 'm2-375EUk', 'embl': '/data/revseq/results/gather_results/RVSeqPlate3-m2/sample_m2-375EUk/m2-375EUk.embl.gz'}, {'name': 'm2-3edUBd', 'embl': '/data/revseq/results/gather_results/RVSeqPlate3-m2/sample_m2-3edUBd/m2-3edUBd.embl.gz'}, {'name': 'm2-kwbwf8', 'embl': '/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-kwbwf8/m2-kwbwf8.embl.gz'}, {'name': 'm2-mc3GSV', 'embl': '/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-mc3GSV/m2-mc3GSV.embl.gz'}, {'name': 'm2-XEuJCy', 'embl': '/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-XEuJCy/m2-XEuJCy.embl.gz'}, {'name': 'm2-UqhsDL', 'embl': '/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-UqhsDL/m2-UqhsDL.embl.gz'}, {'name': 'm2-4UF9K6', 'embl': '/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-4UF9K6/m2-4UF9K6.embl.gz'}, {'name': 'm2-CLwScw', 'embl': '/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-CLwScw/m2-CLwScw.embl.gz'}, {'name': 'm2-jtjCjd', 'embl': '/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-jtjCjd/m2-jtjCjd.embl.gz'}, {'name': 'm2-xwkVbK', 'embl': '/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-xwkVbK/m2-xwkVbK.embl.gz'}, {'name': 'm2-ZAyjrq', 'embl': '/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-ZAyjrq/m2-ZAyjrq.embl.gz'}, {'name': 'm2-RyXauM', 'embl': '/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-RyXauM/m2-RyXauM.embl.gz'}, {'name': 'm2-itJLqU', 'embl': '/data/revseq/results/gather_results/RVSeqPlate4-m2/sample_m2-itJLqU/m2-itJLqU.embl.gz'}, {'name': 'm2-xHGR9u', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-xHGR9u/m2-xHGR9u.embl.gz'}, {'name': 'm2-xrRPFj', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-xrRPFj/m2-xrRPFj.embl.gz'}, {'name': 'm2-xxAxC9', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-xxAxC9/m2-xxAxC9.embl.gz'}, {'name': 'm2-qtZxkz', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-qtZxkz/m2-qtZxkz.embl.gz'}, {'name': 'm2-sxkK9R', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-sxkK9R/m2-sxkK9R.embl.gz'}, {'name': 'm2-SKJgwA', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-SKJgwA/m2-SKJgwA.embl.gz'}, {'name': 'm2-ss4Fh8', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-ss4Fh8/m2-ss4Fh8.embl.gz'}, {'name': 'm2-URVR7h', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-URVR7h/m2-URVR7h.embl.gz'}, {'name': 'm2-Sm9GYA', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-Sm9GYA/m2-Sm9GYA.embl.gz'}, {'name': 'm2-xYnaUy', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-xYnaUy/m2-xYnaUy.embl.gz'}, {'name': 'm2-m6243p', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-m6243p/m2-m6243p.embl.gz'}, {'name': 'm2-opfXjQ', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-opfXjQ/m2-opfXjQ.embl.gz'}, {'name': 'm2-q2oXnd', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-q2oXnd/m2-q2oXnd.embl.gz'}, {'name': 'm2-Qb4JVo', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-Qb4JVo/m2-Qb4JVo.embl.gz'}, {'name': 'm2-g33df2', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-g33df2/m2-g33df2.embl.gz'}, {'name': 'm2-kBgpZ4', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-kBgpZ4/m2-kBgpZ4.embl.gz'}, {'name': 'm2-V5iihT', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-V5iihT/m2-V5iihT.embl.gz'}, {'name': 'm2-SQZYHn', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-SQZYHn/m2-SQZYHn.embl.gz'}, {'name': 'm2-MHWQwo', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-MHWQwo/m2-MHWQwo.embl.gz'}, {'name': 'm2-FHPnrx', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-FHPnrx/m2-FHPnrx.embl.gz'}, {'name': 'm2-VA9Yi5', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-VA9Yi5/m2-VA9Yi5.embl.gz'}, {'name': 'm2-PGaVE8', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-PGaVE8/m2-PGaVE8.embl.gz'}, {'name': 'm2-x6NRYc', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-x6NRYc/m2-x6NRYc.embl.gz'}, {'name': 'm2-oLXZ9a', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-oLXZ9a/m2-oLXZ9a.embl.gz'}, {'name': 'm2-ydSpNQ', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-ydSpNQ/m2-ydSpNQ.embl.gz'}, {'name': 'm2-DwiwwU', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-DwiwwU/m2-DwiwwU.embl.gz'}, {'name': 'm2-39kc9U', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-39kc9U/m2-39kc9U.embl.gz'}, {'name': 'm2-4rvG7f', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-4rvG7f/m2-4rvG7f.embl.gz'}, {'name': 'm2-9HhvMh', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-9HhvMh/m2-9HhvMh.embl.gz'}, {'name': 'm2-eAr6Lo', 'embl': '/data/revseq/results/gather_results/RVSeqPlate5-m2/sample_m2-eAr6Lo/m2-eAr6Lo.embl.gz'}, {'name': 'm2-3VC6ft', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-3VC6ft/m2-3VC6ft.embl.gz'}, {'name': 'm2-sr4bvM', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-sr4bvM/m2-sr4bvM.embl.gz'}, {'name': 'm2-fhoncL', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-fhoncL/m2-fhoncL.embl.gz'}, {'name': 'm2-BnbHLr', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-BnbHLr/m2-BnbHLr.embl.gz'}, {'name': 'm2-fafZrX', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-fafZrX/m2-fafZrX.embl.gz'}, {'name': 'm2-hzNQRG', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-hzNQRG/m2-hzNQRG.embl.gz'}, {'name': 'm2-Rtjhba', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-Rtjhba/m2-Rtjhba.embl.gz'}, {'name': 'm2-LCBdav', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-LCBdav/m2-LCBdav.embl.gz'}, {'name': 'm2-K2udH5', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-K2udH5/m2-K2udH5.embl.gz'}, {'name': 'm2-SpXNmN', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-SpXNmN/m2-SpXNmN.embl.gz'}, {'name': 'm2-4t7Ncr', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-4t7Ncr/m2-4t7Ncr.embl.gz'}, {'name': 'm2-LosGHn', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-LosGHn/m2-LosGHn.embl.gz'}, {'name': 'm2-iDqyV9', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-iDqyV9/m2-iDqyV9.embl.gz'}, {'name': 'm2-voFSf6', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-voFSf6/m2-voFSf6.embl.gz'}, {'name': 'm2-VT3bcF', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-VT3bcF/m2-VT3bcF.embl.gz'}, {'name': 'm2-juVwzj', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-juVwzj/m2-juVwzj.embl.gz'}, {'name': 'm2-gkcF5C', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-gkcF5C/m2-gkcF5C.embl.gz'}, {'name': 'm2-mmYEfJ', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-mmYEfJ/m2-mmYEfJ.embl.gz'}, {'name': 'm2-h2ZJd5', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-h2ZJd5/m2-h2ZJd5.embl.gz'}, {'name': 'm2-uTpHe5', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-uTpHe5/m2-uTpHe5.embl.gz'}, {'name': 'm2-3H4uLZ', 'embl': '/data/revseq/results/gather_results/RVSeqPlate6-m2/sample_m2-3H4uLZ/m2-3H4uLZ.embl.gz'}, {'name': 'm2-Y4wzR3', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-Y4wzR3/m2-Y4wzR3.embl.gz'}, {'name': 'm2-PZQaiq', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-PZQaiq/m2-PZQaiq.embl.gz'}, {'name': 'm2-gF3E5u', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-gF3E5u/m2-gF3E5u.embl.gz'}, {'name': 'm2-tDscqy', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-tDscqy/m2-tDscqy.embl.gz'}, {'name': 'm2-C983Zw', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-C983Zw/m2-C983Zw.embl.gz'}, {'name': 'm2-YBrMoY', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-YBrMoY/m2-YBrMoY.embl.gz'}, {'name': 'm2-gsuEGK', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-gsuEGK/m2-gsuEGK.embl.gz'}, {'name': 'm2-gyUNUP', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-gyUNUP/m2-gyUNUP.embl.gz'}, {'name': 'm2-hPpCR3', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-hPpCR3/m2-hPpCR3.embl.gz'}, {'name': 'm2-HtiJwT', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-HtiJwT/m2-HtiJwT.embl.gz'}, {'name': 'm2-kNU7js', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-kNU7js/m2-kNU7js.embl.gz'}, {'name': 'm2-phtHeq', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-phtHeq/m2-phtHeq.embl.gz'}, {'name': 'm2-jHbrsp', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-jHbrsp/m2-jHbrsp.embl.gz'}, {'name': 'm2-JFKsEJ', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-JFKsEJ/m2-JFKsEJ.embl.gz'}, {'name': 'm2-QqAEGg', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-QqAEGg/m2-QqAEGg.embl.gz'}, {'name': 'm2-XMdsys', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-XMdsys/m2-XMdsys.embl.gz'}, {'name': 'm2-YowhiU', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-YowhiU/m2-YowhiU.embl.gz'}, {'name': 'm2-YYRzxL', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-YYRzxL/m2-YYRzxL.embl.gz'}, {'name': 'm2-hozdkG', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-hozdkG/m2-hozdkG.embl.gz'}, {'name': 'm2-2b8YpQ', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-2b8YpQ/m2-2b8YpQ.embl.gz'}, {'name': 'm2-DgRuyg', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-DgRuyg/m2-DgRuyg.embl.gz'}, {'name': 'm2-PGoRgN', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-PGoRgN/m2-PGoRgN.embl.gz'}, {'name': 'm2-AcEG8S', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-AcEG8S/m2-AcEG8S.embl.gz'}, {'name': 'm2-DMLKnU', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-DMLKnU/m2-DMLKnU.embl.gz'}, {'name': 'm2-os283P', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-os283P/m2-os283P.embl.gz'}, {'name': 'm2-NESNdH', 'embl': '/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-NESNdH/m2-NESNdH.embl.gz'}, {'name': 'm2-b99ooj', 'embl': '/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-b99ooj/m2-b99ooj.embl.gz'}, {'name': 'm2-j3ZXpL', 'embl': '/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-j3ZXpL/m2-j3ZXpL.embl.gz'}, {'name': 'm2-FXFt6J', 'embl': '/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-FXFt6J/m2-FXFt6J.embl.gz'}, {'name': 'm2-LYceeY', 'embl': '/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-LYceeY/m2-LYceeY.embl.gz'}, {'name': 'm2-WrGEQp', 'embl': '/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-WrGEQp/m2-WrGEQp.embl.gz'}, {'name': 'm2-TPsF65', 'embl': '/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-TPsF65/m2-TPsF65.embl.gz'}, {'name': 'm2-f8rfBe', 'embl': '/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-f8rfBe/m2-f8rfBe.embl.gz'}, {'name': 'm2-yhE88h', 'embl': '/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-yhE88h/m2-yhE88h.embl.gz'}, {'name': 'm2-gaQkyT', 'embl': '/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-gaQkyT/m2-gaQkyT.embl.gz'}, {'name': 'm2-FPVUGC', 'embl': '/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-FPVUGC/m2-FPVUGC.embl.gz'}, {'name': 'm2-uSTYZ2', 'embl': '/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-uSTYZ2/m2-uSTYZ2.embl.gz'}, {'name': 'm2-sAyJww', 'embl': '/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-sAyJww/m2-sAyJww.embl.gz'}, {'name': 'm2-2hKeW2', 'embl': '/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-2hKeW2/m2-2hKeW2.embl.gz'}, {'name': 'm2-4Q8Stt', 'embl': '/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-4Q8Stt/m2-4Q8Stt.embl.gz'}, {'name': 'm2-BbET5x', 'embl': '/data/revseq/results/gather_results/RVSeqPlate8-m2/sample_m2-BbET5x/m2-BbET5x.embl.gz'}, {'name': 'm2-6z2wf2', 'embl': '/data/revseq/results/gather_results/RVSeqPlate9-m2/sample_m2-6z2wf2/m2-6z2wf2.embl.gz'}, {'name': 'm2-5cJWyf', 'embl': '/data/revseq/results/gather_results/RVSeqPlate9-m2/sample_m2-5cJWyf/m2-5cJWyf.embl.gz'}, {'name': 'm2-sFFUfJ', 'embl': '/data/revseq/results/gather_results/RVSeqPlate9-m2/sample_m2-sFFUfJ/m2-sFFUfJ.embl.gz'}, {'name': 'm2-TPgies', 'embl': '/data/revseq/results/gather_results/RVSeqPlate9-m2/sample_m2-TPgies/m2-TPgies.embl.gz'}, {'name': 'm2-U2wv8C', 'embl': '/data/revseq/results/gather_results/RVSeqPlate9-m2/sample_m2-U2wv8C/m2-U2wv8C.embl.gz'}, {'name': 'm2-ztFGoy', 'embl': '/data/revseq/results/gather_results/RVSeqPlate10-m2/sample_m2-ztFGoy/m2-ztFGoy.embl.gz'}, {'name': 'm2-3p7s5a', 'embl': '/data/revseq/results/gather_results/RVSeqPlate10-m2/sample_m2-3p7s5a/m2-3p7s5a.embl.gz'}, {'name': 'm2-qeHUBn', 'embl': '/data/revseq/results/gather_results/RVSeqPlate10-m2/sample_m2-qeHUBn/m2-qeHUBn.embl.gz'}, {'name': 'm2-wVBQbT', 'embl': '/data/revseq/results/gather_results/RVSeqPlate11-m2/sample_m2-wVBQbT/m2-wVBQbT.embl.gz'}, {'name': 'm2-bpDXvz', 'embl': '/data/revseq/results/gather_results/RVSeqPlate11-m2/sample_m2-bpDXvz/m2-bpDXvz.embl.gz'}, {'name': 'm2-8ochGa', 'embl': '/data/revseq/results/gather_results/RVSeqPlate13-m2/sample_m2-8ochGa/m2-8ochGa.embl.gz'}, {'name': 'm2-qWgDtq', 'embl': '/data/revseq/results/gather_results/RVSeqPlate13-m2/sample_m2-qWgDtq/m2-qWgDtq.embl.gz'}, {'name': 'm2-33tHNX', 'embl': '/data/revseq/results/gather_results/RVSeqPlate13-m2/sample_m2-33tHNX/m2-33tHNX.embl.gz'}, {'name': 'm2-AcEG8S', 'embl': '/data/revseq/results/gather_results/RVSeqPlate7-m2/sample_m2-AcEG8S/m2-AcEG8S.embl.gz'}, {'name': 'm2-bpDXvz', 'embl': '/data/revseq/results/gather_results/RVSeqPlate11-m2/sample_m2-bpDXvz/m2-bpDXvz.embl.gz'}]


OLD_PREFIX = "/data/revseq/results/gather_results/"
NEW_PREFIX = "/data/"


class Command(BaseCommand):
    help = (
        "Import hardcoded EMBL files stored in the 'data' list above, "
        "normalize the paths, and create File objects linked to Samples."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Only print what would be created, do not touch the DB.",
        )

    def handle(self, *args, **options):
        dry_run = options["dry_run"]

        self.stdout.write(f"Loaded {len(data)} EMBL file entries.")

        # Ensure FileType exists
        embl_type, _ = FileType.objects.get_or_create(postfix="embl")

        # Helper to compute checksum
        def compute_checksum(path_str: str) -> str:
            p = Path(path_str)
            if not p.exists():
                raise FileNotFoundError(f"File does not exist: {path_str}")

            h = hashlib.sha256()
            with p.open("rb") as f:
                for chunk in iter(lambda: f.read(8192), b""):
                    h.update(chunk)
            return h.hexdigest()

        created = 0
        skipped = 0
        errors = 0

        # Dry-run mode (no database transaction)
        if dry_run:
            for item in data:
                name = item["name"]
                raw_path = item["embl"]

                # Normalize path
                normalized = (
                    raw_path.replace(OLD_PREFIX, NEW_PREFIX, 1)
                    if raw_path.startswith(OLD_PREFIX) else raw_path
                )

                # Find sample by pseudonymized_id OR sample_number
                sample = (
                    Sample.objects.filter(pseudonymized_id=name).first()
                    or Sample.objects.filter(sample_number=name).first()
                )
                if not sample:
                    self.stdout.write(self.style.WARNING(
                        f"[DRY RUN] Sample not found for name '{name}'."
                    ))
                    skipped += 1
                    continue

                # If file exists already, skip
                if File.objects.filter(path=normalized).exists():
                    self.stdout.write(
                        f"[DRY RUN] File already exists: {normalized}"
                    )
                    skipped += 1
                    continue

                # Compute checksum
                try:
                    checksum = compute_checksum(normalized)
                except FileNotFoundError as e:
                    self.stdout.write(self.style.ERROR(f"[DRY RUN] {e}"))
                    errors += 1
                    continue

                self.stdout.write(
                    f"[DRY RUN] Would create File:\n"
                    f"          sample={sample.id} ({sample.pseudonymized_id})\n"
                    f"          plate={sample.plate.barcode if sample.plate else 'None'}\n"
                    f"          path={normalized}\n"
                    f"          checksum={checksum[:12]}...\n"
                )

                created += 1

            self.stdout.write(
                self.style.WARNING(
                    f"DRY RUN COMPLETE — would create {created}, skipped {skipped}, errors {errors}."
                )
            )
            return

        # ---------------------------
        # REAL RUN WITH TRANSACTION
        # ---------------------------
        with transaction.atomic():
            for item in data:
                name = item["name"]
                raw_path = item["embl"]

                normalized = (
                    raw_path.replace(OLD_PREFIX, NEW_PREFIX, 1)
                    if raw_path.startswith(OLD_PREFIX) else raw_path
                )

                sample = (
                    Sample.objects.filter(pseudonymized_id=name).first()
                    or Sample.objects.filter(sample_number=name).first()
                )
                if not sample:
                    self.stdout.write(self.style.WARNING(
                        f"Sample not found for name '{name}'. Skipping."
                    ))
                    skipped += 1
                    continue

                if File.objects.filter(path=normalized).exists():
                    self.stdout.write(
                        f"File already exists, skipping: {normalized}"
                    )
                    skipped += 1
                    continue

                try:
                    checksum = compute_checksum(normalized)
                except FileNotFoundError as e:
                    self.stdout.write(self.style.ERROR(str(e)))
                    errors += 1
                    continue

                f = File.objects.create(
                    path=normalized,
                    checksum=checksum,
                    type=embl_type,
                    sample=sample,
                    plate=sample.plate,
                )

                created += 1
                self.stdout.write(self.style.SUCCESS(
                    f"Created File id={f.id} for sample={sample.id} with path={normalized}"
                ))

        self.stdout.write(self.style.SUCCESS(
            f"Done: created {created}, skipped {skipped}, errors {errors}."
        ))