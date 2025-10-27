from django.core.management.base import BaseCommand
from core.models import Sample

"""
This command lists the taxon IDs and scientific names of Substrains associated with a predefined list of pseudonymized sample IDs.

"""
class Command(BaseCommand):
    help = "Prints: sample_id -> Substrain scientific_name (taxon_id) for selected pseudonymized samples."

    # your list of pseudonymized IDs
    SAMPLE_IDS = [
        "ahT4kF", "Q8Ba9A", "33eCP4", "zGrf4n", "HAAUgG", "xPAa2d", "HXGr8V", "wF4k2T",
        "Di2qm7", "iZ7TLc", "3fLBSP", "Kr6366", "5fSMj3", "MT4r8p", "v3KDUF", "Qf8cLC",
        "qFQ9kJ", "3fGJE8", "4JSa55", "G8vc2N", "k9vgiw", "XisgNp", "afGSAT", "CgiSAY",
        "2MHu4K", "TCsq7r", "8YYVHP", "hHKDyB", "G86emu", "6RPWyT", "M5iqYC", "nf9ieP",
        "3eGjGj", "ND22rx", "4i2tZF", "Zu3u7N", "DBKZHX", "FcN3wC", "Ki8Aa3", "4SRYNi",
        "zxhtR4", "T5VXtD", "8ERNNJ", "pzbDPq", "SrvHnQ", "RCdbJH", "sVxWai", "yebiMQ",
        "GbMzJy", "AnXy9n",
    ]

    def handle(self, *args, **options):
        qs = (
            Sample.objects.filter(pseudonymized_id__in=self.SAMPLE_IDS)
            .prefetch_related("samplecounts__substrain")
            .order_by("pseudonymized_id")
        )

        for sample in qs:
            sid = sample.pseudonymized_id
            substrains = []
            seen = set()

            for sc in sample.samplecounts.all():
                ss = sc.substrain
                if ss and ss.taxon_id not in seen:
                    seen.add(ss.taxon_id)
                    name = ss.scientific_name or ss.name or "N/A"
                    substrains.append(f"{name} ({ss.taxon_id or 'N/A'})")

            if substrains:
                print(f"{sid} -> " + "; ".join(substrains))
            else:
                print(f"{sid} -> (no Substrain found)")