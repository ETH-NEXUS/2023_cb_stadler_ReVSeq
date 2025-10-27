from django.core.management.base import BaseCommand
from core.models import Sample

"""
This command lists the taxon IDs and scientific names of Substrains
for a predefined list of pseudonymized sample IDs, sorted by rpkm_proportions (desc),
treating None as -inf.
"""

class Command(BaseCommand):
    help = "Prints: sample_id -> Substrain scientific_name (taxon_id) sorted by rpkm_proportions desc."

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

    @staticmethod
    def _sort_key(sc):
        # None → -inf, so those go to the end when sorting desc
        return sc.rpkm_proportions if sc.rpkm_proportions is not None else float('-inf')

    def handle(self, *args, **options):
        qs = (
            Sample.objects.filter(pseudonymized_id__in=self.SAMPLE_IDS)
            .prefetch_related("samplecounts__substrain")
            .order_by("pseudonymized_id")
        )

        for sample in qs:
            sid = sample.pseudonymized_id
            seen_by_taxon = set()
            items = []

            # sort counts by rpkm_proportions desc (None treated as -inf)
            counts = sorted(sample.samplecounts.all(), key=self._sort_key, reverse=True)

            for sc in counts:
                ss = sc.substrain
                if not ss:
                    continue
                taxon = ss.taxon_id
                if taxon in seen_by_taxon:
                    continue
                seen_by_taxon.add(taxon)
                name = ss.scientific_name or ss.name or "N/A"
                items.append(f"{name} ({taxon if taxon is not None else 'N/A'})")

            if items:
                print(f"{sid} -> " + "\n".join(items[:2]))
                print('----------------------------------------------------------------------')  # blank line for readability
            else:
                print(f"{sid} -> (no Substrain found)")